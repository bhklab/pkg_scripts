library(data.table)
library(CoreGx)
library(PharmacoGx)

if (!require("synergyfinder", quiet = TRUE)) {
    BiocManager::install("synergyfinder")
    library(synergyfinder)
} else {
    library(synergyfinder)
}

if (!require("corrplot", quiet = TRUE)) {
    install.packages("corrplot")
    library(corrplot)
} else {
    library(corrplot)
}

## == Set up parallisation ====================================================
nthread <- 6
setDTthreads(nthread)
## ============================================================================

## == Download the Mathews Griner Drug Combination Screening Dataset ==========
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
## ============================================================================

## == Prepare SynergyFinder input =============================================
## average viability over replicates
viability |> aggregate(
    viability = mean(Value),
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
setkeyv(synergyfinder_input, c("drug_row", "drug_col", "conc_r", "conc_c"))
## ============================================================================

## == Run SynergyFinder =======================================================
bench::system_time({
synergyfinder_input_reshape <- ReshapeData(
    data = as.data.frame(synergyfinder_input),
    data_type = "viability",
    impute = TRUE,
    impute_method = NULL,
    noise = FALSE
)
})
bench::system_time({
synergyfinder_result <- CalculateSynergy(
    data = synergyfinder_input_reshape,
    method = c("ZIP", "HSA", "Bliss", "Loewe"),
    Emin = NA,
    Emax = NA,
    correct_baseline = "non"
)
})
## ============================================================================

## == Run SynergyFinder with Baseline correction & Double Concentrations ======
#bench::system_time({
#synergyfinder_input_reshape_2 <- ReshapeData(
#    data = as.data.frame(copy(synergyfinder_input)[, `:=`(
#        conc_r = conc_r * 2, conc_c = conc_c * 2
#    )]),
#    data_type = "viability",
#    impute = TRUE,
#    impute_method = NULL,
#    noise = FALSE
#)})
#bench::system_time({
#synergyfinder_result_2 <- CalculateSynergy(
#    data = synergyfinder_input_reshape_2,
#    method = c("ZIP", "HSA", "Bliss", "Loewe"),
#    Emin = NA,
#    Emax = NA,
#    correct_baseline = "all"
#)
#})
## ============================================================================

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
pgx_input[, tech_rep := seq_len(.N),
    by = .(treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid)
]

mono_viability_1 <- pgx_input[treatment2dose == 0][, `:=`(
    treatment2id = NULL,
    treatment2dose = NULL
)]

## We consider repeated single drug screening as technical replicates
## i.e. treatment2dose == 0 in each drug combination screening
mono_viability_1 |> aggregate(
    viability = mean(viability),
    by = c("treatment1id", "treatment1dose", "sampleid")
) -> mono_viability_1

mono_viability_1 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment1dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = TRUE),
    enlist = FALSE,
    by = c("treatment1id", "sampleid"),
    nthread = nthread
) -> mono_profiles_1

mono_viability_2 <- pgx_input[treatment1dose == 0][, `:=`(
    treatment1id = NULL,
    treatment1dose = NULL
)]

## We consider repeated single drug screening as technical replicates
## i.e. treatment1dose == 0 in each drug combination screening
mono_viability_2 |> aggregate(
    viability = mean(viability),
    by = c("treatment2id", "treatment2dose", "sampleid")
) -> mono_viability_2

mono_viability_2 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment2dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = TRUE),
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

combo_viability <- pgx_input[
    treatment1dose > 0 & treatment2dose > 0
][, `:=`(
    combo_viability = viability,
    viability = NULL
)]

combo_profiles <- combo_viability[
    mono_profiles_1,
    on = c("treatment1id" = "treatment1id",
           "treatment1dose" = "treatment1dose",
           "sampleid" = "sampleid")
]

combo_profiles <- merge.data.table(
    x = combo_profiles,
    y = mono_profiles_2,
    by = c("treatment2id", "treatment2dose", "sampleid"),
    suffixes = c("_1", "_2")
)
setkeyv(combo_profiles,
        c("treatment1id", "treatment1dose", "tech_rep",
          "treatment2id", "treatment2dose", "sampleid"))
## ============================================================================

## == Calculate Synergy Scores using PharmacoGx ===============================
combo_profiles |> aggregate(
        HSA_ref = 100 - computeHSA(viability_1, viability_2),
        Bliss_ref = 100 - computeBliss(viability_1 / 100, viability_2 / 100) * 100,
        ZIP_ref = 100 - computeZIP(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1 / 100,
            E_inf_2 = E_inf_2 / 100,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        ) * 100,
        Loewe_ref = 100 - PharmacoGx::computeLoewe(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1 / 100,
            E_inf_2 = E_inf_2 / 100,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        ) * 100,
        combo_response = 100 - combo_viability,
        by = key(combo_profiles)
    ) -> pgx_synergy_refs
setkeyv(pgx_synergy_refs,
        c("treatment1id", "treatment1dose",
          "treatment2id", "treatment2dose", "tech_rep"))

pgx_synergy_refs |> aggregate(
    HSA_synergy = combo_response - HSA_ref,
    Bliss_synergy = combo_response - Bliss_ref,
    Loewe_synergy = combo_response - Loewe_ref,
    #HSA_CI = combo_viability / HSA,
    #Bliss_CI = combo_viability / Bliss,
    #Loewe_CI = loeweCI(
    #    viability = combo_viability,
    #    treatment1dose = treatment1dose,
    #    treatment2dose = treatment2dose,
    #    HS_1 = HS_1,
    #    E_inf_1 = E_inf_1,
    #    EC50_1 = EC50_1,
    #    HS_2 = HS_2,
    #    E_inf_2 = E_inf_2,
    #    EC50_2 = EC50_2),
    by = key(combo_profiles),
    nthread = nthread
) -> pgx_synergy_scores

combo_profiles_rep <- copy(combo_profiles)[tech_rep > 1,
    treatment1id := paste0(treatment1id, "_1")]
combo_profiles_rep[,
    .computeZIPdelta(
        treatment1id = treatment1id,
        treatment2id = treatment2id,
        treatment1dose = treatment1dose,
        treatment2dose = treatment2dose,
        sampleid = sampleid,
        HS_1 = HS_1, HS_2 = HS_2,
        EC50_1 = EC50_1, EC50_2 = EC50_2,
        E_inf_1 = E_inf_1 / 100, E_inf_2 = E_inf_2 / 100,
        combo_viability = combo_viability / 100,
        residual = "logcosh",
        #residual = "normal",
        #residual = "Cauchy",
        nthread = nthread,
        show_Rsqr = TRUE
)] -> pgx_delta_scores
pgx_delta_scores[, `:=`(ZIP_synergy = -1 * delta_score * 100, delta_score = NULL)]
pgx_delta_scores[
    grepl(".*_1", treatment1id),
    treatment1id := gsub("_1", "", treatment1id)
]
pgx_delta_scores[, tech_rep := seq_len(.N),
    by = .(treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid)
]

pgx_synergy_scores <- pgx_synergy_scores[pgx_delta_scores, on = c(
    "treatment1id" = "treatment1id",
    "treatment1dose" = "treatment1dose",
    "treatment2id" = "treatment2id",
    "treatment2dose" = "treatment2dose",
    "sampleid" = "sampleid",
    "tech_rep" = "tech_rep"
)][, `:=`(delta_Rsqr_1_to_2 = NULL, delta_Rsqr_2_to_1 = NULL)]

setkeyv(pgx_synergy_scores,
        c("treatment1id", "treatment1dose",
          "treatment2id", "treatment2dose", "tech_rep"))
## ============================================================================

## == PharmacoGx vs. SynergyFinder ============================================
synergyfinder_drug_pairs <- as.data.table(synergyfinder_result$drug_pairs)
synergyfinder_response <- as.data.table(synergyfinder_result$response)
synergyfinder_synergy_scores <- as.data.table(synergyfinder_result$synergy_scores)

synergyfinder_synergy_scores <- synergyfinder_synergy_scores[conc1 > 0 & conc2 > 0]

synergyfinder_synergy_scores <- synergyfinder_synergy_scores[
    synergyfinder_drug_pairs[, .(block_id, drug1, drug2)],
    on = c("block_id" = "block_id")
][, block_id := NULL][, tech_rep := seq_len(.N),
    by = .(drug1, drug2, conc1, conc2)
]
setkeyv(synergyfinder_synergy_scores, c("drug1", "conc1", "drug2", "conc2", "tech_rep"))

synergyfinder_synergy_refs <- synergyfinder_synergy_scores[,
    .(drug1, drug2, conc1, conc2, ZIP_ref, HSA_ref, Bliss_ref, Loewe_ref)
][, tech_rep := seq_len(.N),
    by = .(drug1, drug2, conc1, conc2)
]
setkeyv(synergyfinder_synergy_refs, c("drug1", "conc1", "drug2", "conc2", "tech_rep"))

synergy_refs <- merge.data.table(
    x = pgx_synergy_refs,
    y = synergyfinder_synergy_refs,
    by.x = key(pgx_synergy_refs),
    by.y = key(synergyfinder_synergy_refs),
    suffixes = c("_pgx", "_sf")
)[, combo_response := NULL]

synergy_scores <- merge.data.table(
    x = pgx_synergy_scores[, c(key(pgx_synergy_scores),
        "HSA_synergy", "Bliss_synergy", "Loewe_synergy", "ZIP_synergy"
    ), with = FALSE],
    y = synergyfinder_synergy_scores[, c(key(synergyfinder_synergy_scores),
        "HSA_synergy", "Bliss_synergy", "Loewe_synergy", "ZIP_synergy"
    ), with = FALSE],
    by.x = key(pgx_synergy_scores),
    by.y = key(synergyfinder_synergy_scores),
    suffixes = c("_pgx", "_sf")
)
## ============================================================================

## == Predicted Response Correlation ==========================================
cor_refs <- cor(synergy_refs[, .(
    HSA_ref_pgx, Bliss_ref_pgx, Loewe_ref_pgx, ZIP_ref_pgx,
    HSA_ref_sf, Bliss_ref_sf, Loewe_ref_sf, ZIP_ref_sf
)], method = "pearson", use = "complete.obs")

corrplot(cor_refs, method = "number", type = "upper")
## ============================================================================

## == Synergy Score Correlation ===============================================
cor_scores <- cor(synergy_scores[, .(
    HSA_synergy_pgx, Bliss_synergy_pgx, Loewe_synergy_pgx, ZIP_synergy_pgx,
    HSA_synergy_sf, Bliss_synergy_sf, Loewe_synergy_sf, ZIP_synergy_sf
)], method = "pearson", use = "complete.obs")

corrplot(cor_scores, method = "number", type = "upper")
## ============================================================================
