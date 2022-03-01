library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)

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
Sys.setenv("VROOM_CONNECTION_SIZE" = 10^7)

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
