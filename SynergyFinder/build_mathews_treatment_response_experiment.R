library(data.table)
library(CoreGx)

if (!require("synergyfinder", quiet = TRUE)) {
    BiocManager::install("synergyfinder")
} else {
    library(synergyfinder)
}

## == Set up parallisation ====================================================
nthread <- 6
setDTthreads(nthread)
## ============================================================================

## Data directory
if (!dir.exists("../raw_data"))
    dir.create(path = "../raw_data")

## Download dataset
if (!file.exists("../raw_data/Griner_et_al.zip")) {
    utils::download.file(
        url = "https://ars.els-cdn.com/content/image/1-s2.0-S2001037015000422-mmc4.zip",
        method = "wget",
        destfile = "../raw_data/Griner_et_al.zip"
    )
}

if (!all(file.exists(
    c("../raw_data/metadata.csv",
      "../raw_data/responses.csv")
    ))) {
    utils::unzip(
        zipfile = "../raw_data/Griner_et_al.zip",
        files = c("metadata.csv", "responses.csv"), ## extract dataset files only
        exdir = "../raw_data"
    )
}


viability <- fread("../raw_data/responses.csv")
drug_data <- fread("../raw_data/metadata.csv")

## == Prepare SynergyFinder input =============================================
## average viability over replicates
viability |>
    aggregate(
        viability = mean(Value) / 100,
        by = c("BlockId", "Col", "Row")
    ) -> viability

## "melt" dose-response matrix to a long-format table
dt_blocks <- mapply(
    FUN = function(x, y) {
        n_row <- max(x$Row)
        n_col <- max(x$Col)
        dt_block <- data.table(
            block_id = y$BlockId,
            drug_row = rep(y$RowName, each = n_row),
            drug_col = rep(y$ColName, n_col),
            conc_r = rep(as.numeric(unlist(strsplit(y$RowConcs, split = ","))), each = n_row),
            conc_c = rep(as.numeric(unlist(strsplit(y$ColConcs, split = ","))), n_col),
            response = x$viability,
            conc_r_unit = rep(y$RowConcUnit, each = n_row),
            conc_c_unit = rep(y$ColConcUnit, n_col)
        )
        return(dt_block)
    },
    split(viability, by = "BlockId"),
    split(drug_data, by = "BlockId"),
    SIMPLIFY = FALSE
)
synergyfinder_input <- rbindlist(dt_blocks)
## == Run SynergyFinder =======================================================
synergyfinder_input_reshape <- ReshapeData(
    data = as.data.frame(synergyfinder_input),
    data_type = "viability",
    impute = TRUE,
    impute_method = NULL,
    noise = FALSE
)
synergyfinder_result <- CalculateSynergy(
    data = synergyfinder_input_reshape,
    method = c("ZIP", "HSA", "Bliss", "Loewe"),
    Emin = NA,
    Emax = NA,
    correct_baseline = "non"
)
## == Prepare PharmacoGx input ================================================
pgx_input <- copy(synergyfinder_input)
pgx_input[, `:=`(
    block_id = NULL,
    sampleid = "DLBCL",
    conc_r_unit = NULL,
    conc_c_unit = NULL
)]
setnames(
    pgx_input,
    old = c("drug_row", "drug_col", "conc_r", "conc_c", "response"),
    new = c("treatment1id",
            "treatment2id",
            "treatment1dose",
            "treatment2dose",
            "viability")
)

mono_viability_1 <- pgx_input[treatment2dose == 0][, `:=`(
    treatment2id = NULL,
    treatment2dose = NULL
)]

mono_viability_1 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment1dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = FALSE),
    enlist = FALSE,
    by = c("treatment1id", "sampleid"),
    nthread = nthread
) -> mono_profiles_1

mono_viability_2 <- pgx_input[treatment1dose == 0][, `:=`(
    treatment1id = NULL,
    treatment1dose = NULL
)]

mono_viability_2 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment2dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = FALSE),
    enlist = FALSE,
    by = c("treatment2id", "sampleid"),
    nthread = 1
) -> mono_profiles_2

mono_viability_1[mono_profiles_1, on = c(
    "treatment1id" = "treatment1id",
    "sampleid" = "sampleid"
)] -> mono_profiles_1

mono_viability_2[mono_profiles_2, on = c(
    "treatment2id" = "treatment2id",
    "sampleid" = "sampleid"
)] -> mono_profiles_2

