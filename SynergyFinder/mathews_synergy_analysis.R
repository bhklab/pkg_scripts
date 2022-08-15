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
        ## Correct baseline using viability administered with 0 concentration
        baseline <- x$viability[which(x$Col == 6 & x$Row == 6)]
        dt_block <- data.table(
            block_id = y$BlockId,
            drug_row = rep(y$RowName, each = n_row),
            drug_col = rep(y$ColName, n_col),
            conc_r = rep(as.numeric(unlist(strsplit(y$RowConcs, split = ","))), each = n_row),
            conc_c = rep(as.numeric(unlist(strsplit(y$ColConcs, split = ","))), n_col),
            response = (x$viability / baseline) * 100,
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
    #block_id = NULL,
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
pgx_input[, bio_rep := seq_len(.N),
    by = .(treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid)
]

mono_viability_1 <- pgx_input[treatment2dose == 0][, `:=`(
    treatment2id = NULL,
    treatment2dose = NULL
)]

mono_viability_1 |> aggregate(
    viability = mean(viability),
    by = c("block_id", "treatment1id", "treatment1dose", "sampleid", "bio_rep")
) -> mono_viability_1

mono_viability_1 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment1dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = TRUE,
        #family = "Cauchy"
        family = "normal"
        ),
    enlist = FALSE,
    by = c("block_id", "treatment1id", "sampleid", "bio_rep"),
    nthread = nthread
) -> mono_profiles_1

mono_viability_2 <- pgx_input[treatment1dose == 0][, `:=`(
    treatment1id = NULL,
    treatment1dose = NULL
)]

mono_viability_2 |> aggregate(
    viability = mean(viability),
    by = c("block_id", "treatment2id", "treatment2dose", "sampleid", "bio_rep")
) -> mono_viability_2

mono_viability_2 |> aggregate(
    PharmacoGx::logLogisticRegression(
        conc = treatment2dose,
        viability = viability,
        trunc = FALSE,
        viability_as_pct = TRUE,
        family = "normal"
        #family = "Cauchy"
    ),
    enlist = FALSE,
    by = c("block_id", "treatment2id", "sampleid", "bio_rep"),
    nthread = 1
) -> mono_profiles_2

mono_viability_1[mono_profiles_1, on = c(
    "block_id" = "block_id",
    "treatment1id" = "treatment1id",
    "sampleid" = "sampleid",
    "bio_rep" = "bio_rep"
)] -> mono_profiles_1

mono_viability_2[mono_profiles_2, on = c(
    "block_id" = "block_id",
    "treatment2id" = "treatment2id",
    "sampleid" = "sampleid",
    "bio_rep" = "bio_rep"
)] -> mono_profiles_2

combo_viability <- pgx_input[
    treatment1dose > 0 & treatment2dose > 0
][, `:=`(
    combo_viability = viability,
    viability = NULL
)]

combo_profiles <- combo_viability[
    mono_profiles_1,
    on = c("block_id" = "block_id",
           "treatment1id" = "treatment1id",
           "treatment1dose" = "treatment1dose",
           "sampleid" = "sampleid",
           "bio_rep" = "bio_rep")
]

combo_profiles <- merge.data.table(
    x = combo_profiles,
    y = mono_profiles_2,
    by = c("block_id" = "block_id",
           "treatment2id",
           "treatment2dose",
           "sampleid",
           "bio_rep"),
    suffixes = c("_1", "_2")
)
setkeyv(combo_profiles,
        c("treatment1id", "treatment1dose", "bio_rep",
          "treatment2id", "treatment2dose", "sampleid"))
## ============================================================================

## == Calculate Synergy Scores using PharmacoGx ===============================
combo_profiles |> aggregate(
        HSA_ref = 100 - computeHSA(viability_1, viability_2),
        Bliss_ref = (1 - computeBliss(viability_1 / 100, viability_2 / 100)) * 100,
        ZIP_ref = 100 * (1 - computeZIP(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1 / 100,
            E_inf_2 = E_inf_2 / 100,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        )),
        Loewe_ref = 100 * (1 - PharmacoGx::computeLoewe(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1 / 100,
            E_inf_2 = E_inf_2 / 100,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        )),
        combo_response = 100 - combo_viability,
        by = c("block_id", key(combo_profiles))
    ) |> na.omit() -> pgx_synergy_refs
setkeyv(pgx_synergy_refs,
        c("block_id", "treatment1id", "treatment1dose",
          "treatment2id", "treatment2dose", "bio_rep"))

combo_profiles_rep <- copy(combo_profiles)[bio_rep > 1,
    treatment1id := paste0(treatment1id, "_1")][,
        `:=`(
            E_inf_1 = E_inf_1 / 100,
            E_inf_2 = E_inf_2 / 100,
            viability_1 = viability_1 / 100,
            viability_2 = viability_2 / 100,
            combo_viability = combo_viability / 100
        )
    ]

## Fit two-way projected Hill curves
ZIP_fit_pgx <- fitTwowayZIP(combo_profiles_rep, residual = "logcosh")
ZIP_fit_pgx[,
    ZIP_fit := ( 1 - 
        hillCurve(
            dose = log10(treatment1dose),
            EC50 = log10(EC50_proj_2_to_1),
            HS = HS_proj_2_to_1,
            E_inf = E_inf_proj_2_to_1,
            E_ninf = E_ninf_proj_2_to_1
        ) + 1 -
        hillCurve(
            dose = log10(treatment2dose),
            EC50 = log10(EC50_proj_1_to_2),
            HS = HS_proj_1_to_2,
            E_inf = E_inf_proj_1_to_2,
            E_ninf = E_ninf_proj_1_to_2
        )
    ) * 100 / 2
]
ZIP_fit_pgx[, bio_rep := seq_len(.N),
    by = .(treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid)
]
ZIP_fit_pgx[
    grepl(".*_1", treatment1id),
    `:=`(treatment1id = gsub("_1", "", treatment1id), bio_rep = bio_rep + 1)
]

pgx_synergy_refs <- pgx_synergy_refs[
    ZIP_fit_pgx[,c(key(pgx_synergy_refs), "ZIP_fit"), with = FALSE],
    on = c(
        "block_id" = "block_id",
        "treatment1id" = "treatment1id",
        "treatment1dose" = "treatment1dose",
        "treatment2id" = "treatment2id",
        "treatment2dose" = "treatment2dose",
        "bio_rep" = "bio_rep"
    )
] |> na.omit()
setkeyv(pgx_synergy_refs,
        c("block_id", "treatment1id", "treatment1dose",
          "treatment2id", "treatment2dose", "bio_rep"))

pgx_synergy_refs |> aggregate(
    HSA_synergy = combo_response - HSA_ref,
    Bliss_synergy = combo_response - Bliss_ref,
    Loewe_synergy = combo_response - Loewe_ref,
    ZIP_synergy = ZIP_fit - ZIP_ref,
    by = c(key(pgx_synergy_refs)),
    nthread = nthread
) -> pgx_synergy_scores

## Alternatively, one can use the vectorised version of computeZIPdelta
## However, since we are converting from viability to response
## to align with the effect measurement used by SynergyFinder,
## Splitting into multiple steps is a safer approach.
#combo_profiles_rep[,
#    .computeZIPdelta(
#        treatment1id = treatment1id,
#        treatment2id = treatment2id,
#        treatment1dose = treatment1dose,
#        treatment2dose = treatment2dose,
#        sampleid = sampleid,
#        HS_1 = HS_1, HS_2 = HS_2,
#        EC50_1 = EC50_1, EC50_2 = EC50_2,
#        E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
#        combo_viability = combo_viability,
#        residual = "logcosh",
#        nthread = nthread,
#        show_Rsqr = TRUE
#)] -> pgx_delta_scores
#pgx_delta_scores[, `:=`(ZIP_synergy = -1 * delta_score * 100, delta_score = NULL)]
#pgx_delta_scores[, bio_rep := seq_len(.N),
#    by = .(treatment1id, treatment2id, treatment1dose, treatment2dose, sampleid)
#]
#pgx_delta_scores[
#    grepl(".*_1", treatment1id),
#    `:=`(treatment1id = gsub("_1", "", treatment1id), bio_rep = bio_rep + 1)
#]
#
#pgx_synergy_scores <- pgx_synergy_scores[pgx_delta_scores, on = c(
#    "treatment1id" = "treatment1id",
#    "treatment1dose" = "treatment1dose",
#    "treatment2id" = "treatment2id",
#    "treatment2dose" = "treatment2dose",
#    "bio_rep" = "bio_rep"
#)][, `:=`(delta_Rsqr_1_to_2 = NULL, delta_Rsqr_2_to_1 = NULL)]

setkeyv(pgx_synergy_scores,
        c("block_id", "treatment1id", "treatment1dose",
          "treatment2id", "treatment2dose", "bio_rep"))
## ============================================================================

## == PharmacoGx vs. SynergyFinder ============================================
synergyfinder_drug_pairs <- as.data.table(synergyfinder_result$drug_pairs)
synergyfinder_response <- as.data.table(synergyfinder_result$response)
synergyfinder_synergy_scores <- as.data.table(synergyfinder_result$synergy_scores)

synergyfinder_synergy_scores <- synergyfinder_synergy_scores[conc1 > 0 & conc2 > 0]

synergyfinder_synergy_scores <- synergyfinder_synergy_scores[
    synergyfinder_drug_pairs[, .(block_id, drug1, drug2)],
    on = c("block_id" = "block_id")
][, bio_rep := seq_len(.N),
    by = .(drug1, drug2, conc1, conc2)
]

synergyfinder_synergy_refs <- synergyfinder_synergy_scores[,
    .(block_id, drug1, drug2, conc1, conc2, ZIP_ref, ZIP_fit, HSA_ref, Bliss_ref, Loewe_ref)
][, bio_rep := seq_len(.N),
    by = .(drug1, drug2, conc1, conc2)
]
setkeyv(synergyfinder_synergy_refs, c("block_id", "drug1", "conc1", "drug2", "conc2", "bio_rep"))

synergyfinder_synergy_scores <- synergyfinder_synergy_refs[
    synergyfinder_response,
    on = c(
        "block_id" = "block_id",
        "conc1" = "conc1",
        "conc2" = "conc2"
    )
] |> na.omit()

synergyfinder_synergy_scores |> aggregate(
    HSA_synergy = response - HSA_ref,
    Bliss_synergy = response - Bliss_ref,
    Loewe_synergy = response - Loewe_ref,
    ZIP_synergy = ZIP_fit - ZIP_ref,
    by = c("block_id", "drug1", "drug2", "conc1", "conc2", "bio_rep")
) -> synergyfinder_synergy_scores
setkeyv(synergyfinder_synergy_scores, c("block_id", "drug1", "conc1", "drug2", "conc2", "bio_rep"))

synergy_refs <- merge.data.table(
    x = pgx_synergy_refs,
    y = synergyfinder_synergy_refs,
    by.x = key(pgx_synergy_refs),
    by.y = key(synergyfinder_synergy_refs),
    suffixes = c("_pgx", "_sf")
)

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
    HSA_ref_pgx, Bliss_ref_pgx, Loewe_ref_pgx, ZIP_ref_pgx, ZIP_fit_pgx,
    HSA_ref_sf, Bliss_ref_sf, Loewe_ref_sf, ZIP_ref_sf, ZIP_fit_sf
)], method = "pearson", use = "complete.obs")

corrplot(cor_refs[
            c("HSA_ref_pgx", "Bliss_ref_pgx", "Loewe_ref_pgx", "ZIP_ref_pgx"),
            c("HSA_ref_sf",  "Bliss_ref_sf",  "Loewe_ref_sf",  "ZIP_ref_sf")
         ],
         method = "color",
         type = "upper",
         addCoef.col = 'white',
         cl.pos = 'b'
)

corrplot(cor_refs[
            c("ZIP_fit_pgx", "ZIP_fit_sf"),
            c("ZIP_fit_sf", "ZIP_fit_pgx")
         ],
         method = "color",
         type = "upper",
         addCoef.col = 'white',
         cl.pos = 'b'
)
## ============================================================================

## == Synergy Score Correlation ===============================================
cor_scores <- cor(synergy_scores[, .(
    HSA_synergy_pgx, Bliss_synergy_pgx, Loewe_synergy_pgx, ZIP_synergy_pgx,
    HSA_synergy_sf, Bliss_synergy_sf, Loewe_synergy_sf, ZIP_synergy_sf
)], method = "pearson", use = "complete.obs")

corrplot(cor_scores[
            c("HSA_synergy_pgx", "Bliss_synergy_pgx", "Loewe_synergy_pgx", "ZIP_synergy_pgx"),
            c("HSA_synergy_sf", "Bliss_synergy_sf", "Loewe_synergy_sf", "ZIP_synergy_sf")
         ],
         method = "color",
         type = "upper",
         addCoef.col = 'white',
         cl.pos = 'b'
)
## ============================================================================




