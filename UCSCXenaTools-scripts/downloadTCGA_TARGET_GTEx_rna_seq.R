library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)
library(qs)

## NOTE:: This script requires ~ 40 GB of RAM

# -- Find datasets of interest

# Load dataset metadata index
data(XenaData)
xd <- as.data.table(XenaData)

(xena_hosts <- unique(xd$XenaHostNames))
host <- "toilHub"

(tcga_cohorts <- unique(xd[XenaCohorts %ilike% "^TCGA", ]$XenaCohorts))
cohort <- "TCGA TARGET GTEx"

(datatypes <- xd[XenaCohorts == cohort, unique(DataSubtype)])
datatype <- "gene expression RNAseq"

(labels <- xd[XenaCohorts == cohort & DataSubtype == datatype, unique(Label)])
label <- "RSEM tpm"

(datasets <- xd[
    XenaCohorts == cohort & DataSubtype == datatype & Label == label,
    XenaDatasets
])

# -- Download the phenotype data

# Extend download.file timeout to 10 minutes
ops <- options()
on.exit(options(ops))
options(timeout=600)

# Increase size of VROOM file buffer
Sys.setenv("VROOM_CONNECTION_SIZE" = 10^8)

XenaGenerate(subset=XenaCohorts == cohort & Type == "clinicalMatrix") |>
    XenaQuery() |>
    XenaDownload(destdir="local_data", download_probeMap=TRUE) |>
    XenaPrepare() ->
    pheno_data

# -- Download the molecular data
XenaGenerate(subset=XenaCohorts == cohort & XenaDatasets %in% datasets) |>
    XenaQuery() |>
    XenaDownload(destdir="local_data", download_probeMap=TRUE) |>
    XenaPrepare() ->
    mol_data

# -- Build MultiAssayExperiment

# colData
for (df in pheno_data) setDT(df)
pheno_idx <- grep("phenotype", names(pheno_data))
cData <- copy(pheno_data[[pheno_idx]])
setkeyv(cData, "sample")
pData <- copy(pheno_data[-pheno_idx])
for (df in pData) {
    setkeyv(df, "sample")
    cData <- df[cData, , on="sample"]  # data.table subset join x[y, ]
}

# rowData
for (df in mol_data) setDT(df)
probemap_idx <- grep("probemap", names(mol_data))
rData <- copy(mol_data[[probemap_idx]])
mData <- copy(mol_data[-probemap_idx])

# assays
for (i in seq_along(mData)) {
    mData[[i]] <- as.matrix(mData[[i]], rownames="feature")
}

# metadata
metadata <- list(
    phenotype_data=xd[XenaCohorts == cohort & Type == "clinicalMatrix"],
    molecular_data=xd[XenaDatasets %in% datasets, ],
    sessionInfo=sessionInfo()
)

rm(mol_data, gc=TRUE)

# experiments
se_list <- vector("list", length(mData))
setkeyv(cData, "sample")
setkeyv(rData, "id")
for (i in seq_along(mData)) {
    se_list[[i]] <- SummarizedExperiment(
        assay=mData[[i]],
        colData=cData[colnames(mData[[i]]), ],
        rowData=rData[rownames(mData[[i]]), ]
    )
}
se_list <- setNames(se_list, names(mData))

# MultiAssayExperient
mae <- MultiAssayExperiment(experiments=se_list, metadata=metadata)

message("INFO: This object is ", format(object.size(mae), "GiB"), " in memory!")
qsave(mae, file=file.path("local_data", paste0(make.names(cohort), "_mae.qs")),
    nthread=getDTthreads())