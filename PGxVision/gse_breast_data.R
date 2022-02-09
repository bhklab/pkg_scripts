library(qs)
library(data.table)
library(SummarizedExperiment)

# -- load the breast cancer summarized experiment
data_dir <- file.path(".", "local_data")
gse_file <- file.path(data_dir, "GSE20194_SE.qs")
gse_breast <- qread(gse_file)

# -- parse to data table and save to disk
gse_df <- as.data.table(assays(gse_breast)[[1L]], keep.rownames="feature")
gse_df <- melt.data.table(gse_df, id.vars="feature", variable.name="sample",
    value.name="expression")

fwrite(gse_df, file=file.path(data_dir, "gse_breast_df.csv"))