library(GEOquery)
library(Biobase)
library(SummarizedExperiment)
library(data.table)
library(qs)
library(beadarray)
library(BeadArrayUseCases)

# AnnotationDbi package for Illumina Bead Array nuIDs
if (!require(lumiHumanAll.db)) BiocManager::install("lumiHumanALl.db")
library(lumiHumanAll.db)

# Convert the annotations to a data.table
ensembl_features <- as.data.table(as.data.frame(lumiHumanAllENSEMBL))

## This pipeline was adapted from:
## https://bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf

# configure download options for this session
ops <- options()
options(timeout=1e6)
on.exit(options(ops))

# configure which datasets to download
datasets <- c("GSE21257")
data_dir <- "local_data"

# fetch Gencode v33 annotations from BHKLAB-Pachyderm/Annotations
gencode_url <- "https://github.com/BHKLAB-Pachyderm/Annotations/raw/master/Gencode.v33.annotation.RData"
gencode_file <- file.path(data_dir, "gencode_annot.RData")
download.file(gencode_url, destfile=gencode_file)
gencode_names <- load(gencode_file)
gencode_annots <- lapply(gencode_names, get) |> setNames(gencode_names)
gene_annots <- gencode_annots$features_gene
tx_annots <- gencode_annots$features_transcript

# clean up gene_ids to match array
setDT(gene_annots)
gene_annots[,
    c("gene_id_versioned", "gene_id") := .(gene_id, gsub("\\..*$", "", gene_id,))
]
setkeyv(gene_annots, "gene_id")

se_list <- vector("list", length(datasets)) |> setNames(datasets)
for (ds in datasets) {
    print(ds)
    eset <- getGEO(ds)[[1]]
    fdata <- as.data.table(as(featureData(eset), "data.frame"))
    fdata <- merge.data.table(fdata, ensembl_features[!duplicated(probe_id)],
        by.x="ID", by.y="probe_id", all.x=TRUE, sort=FALSE)
    # drop Y chromosome to prevent duplicated probe_ids due to genes shared
    # with X chromosome
    fdata <- merge.data.table(fdata, gene_annots[seqnames != "chrY", ],
        by.x="ensembl_id", by.y="gene_id", all.x=TRUE, sort=FALSE)
    featData <- AnnotatedDataFrame(fdata)
    row.names(featData) <- fdata$ID
    featureData(eset) <- featureData
    rownames(eset) <- featureData(eset)$ID
    # coerce to SummarizedExperiment
    se_type <- if (any(is.na(featureData(eset)$start))) "SummarizedExperiment"
        else "RangedSummarizedExperiment"
    se_list[[ds]] <- as(eset, se_type)
}

# Save datasets to disk
for (ds in names(se_list)) {
    qsave(se_list[[ds]],
        file=file.path(data_dir, paste0(ds, "_", class(se_list[[ds]])[1], ".qs")),
        nthread=getDTthreads()
    )
}