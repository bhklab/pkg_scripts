library(qs)
library(data.table)

# -- Download the tcga data from ORCESTRA
data_dir <- file.path(".", "local_data")
data_url <- "https://zenodo.org/record/5731017/files/TCGA_BRCA.rds?download=1"
tcga_breast <- file.path(data_dir, "TCGA_BRCA.rds")

# increase download timeout, zenodo is slow
opts <- options()
options(timeout=9000)
on.exit(options(opts))

if (!file.exits(tcga_breast)) download.file(url=data_url, destfile=tcga_breast)

# -- Load the dataset, extract microarray data, parse to flat file
tcga_mae <- readRDS(tcga_breast)

assay_name <- "BRCA_mRNAArray-20160128"
tcga_se <- getWithColData(tcga_mae, assay_name)

replicates <- replicated(tcga_mae)[[assay_name]]
tcga_se_no_reps <- mergeReplicates(tcga_se, replicates)

tcga_df <- as.data.table(assays(tcga_se)[[1L]], keep.rownames="feature")
tcga_df <- melt.data.table(tcga_df, id.vars="feature", variable.name="sample",
    value.name="expression")

fwrite(tcga_df, file=file.path(".", "local_data", "tcga_breast_df.csv"))