library(qs)
library(data.table)
library(MultiAssayExperiment)

# -- Download the tcga data from ORCESTRA
data_dir <- file.path(".", "local_data")
data_url <- "https://zenodo.org/record/5731017/files/TCGA_BRCA.rds?download=1"
tcga_breast <- file.path(data_dir, "TCGA_BRCA.rds")

# increase download timeout, zenodo is slow
opts <- options()
options(timeout=9000)
on.exit(options(opts))

if (!file.exists(tcga_breast)) download.file(url=data_url, destfile=tcga_breast)

# -- Load the dataset, extract microarray data, parse to flat file
tcga_mae <- readRDS(tcga_breast)

assay_name <- "BRCA_RNASeq2GeneNorm-20160128"
tcga_se <- getWithColData(tcga_mae, assay_name)

# merge replicates by mean
replicates <- replicated(tcga_mae)[[assay_name]]
tcga_se_no_reps <- mergeReplicates(tcga_se, replicates)

tcga_df <- as.data.table(
    round(log2(assays(tcga_se_no_reps)[[1L]] + 1e-6), 2),
    keep.rownames="feature"
)
tcga_df_melt <- melt.data.table(tcga_df, id.vars="feature", variable.name="sample",
    value.name="expression")

fwrite(tcga_df, file=file.path(".", "local_data", "tcga_breast_rnaseq_mat.csv"))

# get interesting patients based on biomarker of interest
set.seed(19900501)
biomarkers <- c("ERBB2", "EGFR", "MET", "GPRC5A")
for (feat in biomarkers) {
    # get random sample over specified percentile
    sample <- tcga_df_melt[
        feature == feat,
        .SD[expression > quantile(expression, 0.90)]
    ][
        sample(seq_len(.N), 1),
        as.character(sample)
    ][1]
    print(sample)
    patient_df <- tcga_df[, c("feature", sample), with=FALSE]
    fwrite(patient_df, file=file.path("local_data",
        paste0("tcga_breast_", feat, "_", sample, ".csv")))
}