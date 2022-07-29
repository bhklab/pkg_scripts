library(data.table)
library(synergyfinder)
library(PharmacoGx)
library(future) ## parallelsing synergyfinder

## == Set up parallisation ====================================================
nthread <- 6
setDTthreads(nthread)
options(mc.cores = nthread)
options(future.globals.maxSize = 1000 * 1024^2) ## 1Gb object size limit
plan(
    strategy = "multicore", ## use multisession if you are running this on Windows
    workers = nthread
)
## ============================================================================

###############################################################################
use_all_single_dose = FALSE
###############################################################################

data_path <- "~/bhklab/feifei/data" ## add your path to NCI PSet here
nci <- readRDS(file.path(data_path, "NCI_ALMANAC_2017.rds"))

bench::system_time({
treatmentResponse(nci) |>
    endoaggregate(
        assay = "sensitivity",
        target = "mono_viability",
        subset = treatment2id == "",
        viability = mean(viability) / 100,
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) -> tre
})

bench::system_time({
tre |>
    endoaggregate(
        assay = "sensitivity",
        target = "combo_viability",
        subset = (!is.na(treatment2dose)),
        combo_viability = mean(viability) / 100,
        by=c("treatment1id", "treatment1dose", "treatment2id", "treatment2dose", "sampleid")
    ) -> tre
})

combo_keys <- c("treatment1id", "treatment1dose", "treatment2id", "treatment2dose", "sampleid")
combo_viability <- tre$combo_viability

## Build a new combo_viability out of mono_viability with one traetment with dose 0
mono_viability <- tre$mono_viability
drug_cell_combo <- combo_viability[, .(treatment1id, treatment2id, sampleid)] |> unique()
combo_viability_0 <- drug_cell_combo ## Can you guess why the 0?
combo_viability_0[, `:=`(treatment1dose = 0, treatment2dose = 0, combo_viability = 1)]

if (use_all_single_dose) {
    ## use all single drug screening dose
    combo_viability_1 <- drug_cell_combo[mono_viability, , on = c(
        "treatment1id" = "treatment1id",
        "sampleid" = "sampleid"
    ), allow.cartesian = TRUE ## I know what I'm doing
    ] 
    
    ## drug combos with treatment2dose = 0
    combo_viability_1[, `:=`(
        treatment2dose = 0,
        combo_viability = viability,
        viability = NULL
    )]
    
    
    combo_viability_2 <- drug_cell_combo[mono_viability, , on = c(
        "treatment2id" = "treatment1id",
        "sampleid" = "sampleid"
    ), allow.cartesian = TRUE ## I know what I'm doing
    ] 
    
    ## drug combos with treatment1dose = 0
    combo_viability_2[, `:=`(
        treatment2dose = treatment1dose,
        treatment1dose = 0,
        combo_viability = viability,
        viability = NULL
    )]
    
    
    ## Remove drugs not screened as treatment 2 in combinations
    combo_viability_1 <- combo_viability_1[!is.na(treatment2id)]
    ## Remove drugs not screened as treatment 1 in combinations
    combo_viability_2 <- combo_viability_2[!is.na(treatment1id)]
} else {
    ## Use only single drug doses present in combination screening
    combo_viability_1 <- combo_viability[, .(treatment1id, treatment1dose,
                                         treatment2id, treatment2dose,
                                         sampleid)][, treatment2dose := 0]

    combo_viability_1 <- combo_viability_1[mono_viability, on = c(
        "treatment1id" = "treatment1id",
        "treatment1dose" = "treatment1dose",
        "sampleid" = "sampleid"
    )][!is.na(treatment2id)][, `:=`(combo_viability = viability, viability = NULL)]
    
    combo_viability_2 <- combo_viability[, .(treatment1id, treatment1dose,
                                         treatment2id, treatment2dose,
                                         sampleid)][, treatment1dose := 0]
    combo_viability_2 <- combo_viability_2[mono_viability, on = c(
        "treatment2id" = "treatment1id",
        "treatment2dose" = "treatment1dose",
        "sampleid" = "sampleid"
    )][!is.na(treatment1id)][, `:=`(combo_viability = viability, viability = NULL)]
}

setcolorder(combo_viability_0, colnames(combo_viability))
setcolorder(combo_viability_1, colnames(combo_viability))
setcolorder(combo_viability_2, colnames(combo_viability))

synergyfinder_input <- rbindlist(
    list(combo_viability_1, combo_viability_2, combo_viability)
)
synergyfinder_input <- unique(synergyfinder_input) ## duplicate rows from joining mono_viability with combo_viability
block_id_cols <- c("treatment1id", "treatment2id", "sampleid")
synergyfinder_input[, block_id := .GRP, by = block_id_cols]

tre_metadata <- metadata(tre)$experiment_metadata

synergyfinder_input[, `:=`(
    conc_unit1 = tre_metadata$treatment1unit,
    conc_unit2 = tre_metadata$treatment2unit
)]

setnames(synergyfinder_input,
         c("treatment1id", "treatment2id",
           "treatment1dose", "treatment2dose",
           "sampleid", "combo_viability"),
         c("drug1", "drug2", "conc1", "conc2", "cellline", "response"))

synergyfinder_input <- as.data.frame(synergyfinder_input)

# run on a single core for debugging
#bench::system_time({
#    synergyfinder_input_reshape <- ReshapeData(data = synergyfinder_input,impute = T)
#}) -> synergyfinder_reshape_time

bench::system_time({
    parallelised_reshape <- future({
        ReshapeData(data = synergyfinder_input,
                    data_type = "viability",
                    impute = TRUE) ## their assumption for the experiment design leaves me no choice
    })
}) -> synergyfinder_reshape_time

synergyfinder_input_reshape <- future::value(parallelised_reshape)

saveRDS(synergyfinder_input_reshape, file.path(data_path, "synergyfinder_input_reshape.rds"))
saveRDS(synergyfinder_reshape_time, file.path(data_path, "synergyfinder_reshape_time.rds"))

bench::system_time({
    parallelised_synergy <- future({
        CalculateSynergy(
            data = synergyfinder_input_reshape,
            method = c("ZIP", "HSA", "Bliss", "Loewe"),
            Emin = NA,
            Emax = NA,
            correct_baseline = "non"
        )
    })
}) -> synergyfinder_compute_time

synergyfinder_result <- future::value(parallelised_synergy)

saveRDS(synergyfinder_result, file.path(data_path, "synergyfinder_result.rds"))
saveRDS(synergyfinder_compute_time, file.path(data_path, "synergyfinder_compute_time.rds"))

