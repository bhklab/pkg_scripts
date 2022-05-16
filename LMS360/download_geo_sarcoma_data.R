library(GEOquery)
library(Biobase)
library(SummarizedExperiment)
library(data.table)
library(qs)
library(affy)

# configure download options for this session
ops <- options()
options(timeout=1e6)
on.exit(options(ops))
# install appropriate brain array annotations for array platform
brain_array_urls <- function(array, species="hs", annotation="ensg", version="25.0.0") {
    ## FIXME:: make robust to missing arguments
    paste0("http://mbni.org/customcdf/", version, "/", annotation, ".download",
        "/", c("", "", "pd."), array, c("", "", "."), species, c("", "", "."),
        annotation, c("cdf", "probe", ""), "_", version, ".tar.gz")
}
array <- "hgu133a"
brain_array <- brain_array_urls(array)
for (pkg in brain_array) {
    install.packages(pkg, type="src", repos=NULL)
}

cdf <- gsub("_.*$", "", basename(brain_array[1]))
library(cdf, character.only=TRUE)

datasets <- c("GSE21122", "GSE21050")
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

# get most recent Affymetrix probe annotations for the HG-U133A array
# Source: "http://www.affymetrix.com/Auth/analysis/downloads/na35/ivt/HG-U133A_2.na35.annot.csv.zip"
# Note: can't download programatically due to account requirements for Affymetrix website
# Note: these are not used, using gencode instead
# affy_probes <- fread(file.path(data_dir, "HG-U133A_2.na35.annot.csv"),
#     skip=25)  # first 25 lines are header

se_list <- vector("list", length(datasets)) |> setNames(datasets)
for (ds in datasets) {
    print(ds)
    eset <- getGEO(ds)
    # extract sample metadata
    pData <- phenoData(eset[[1]])
    # download raw .CEL files
    file_df <- getGEOSuppFiles(ds,
        baseDir=data_dir,
        filter_regex=".*RAW.tar"
    )
    dataset_dir <- dirname(rownames(file_df)[1])
    untar(rownames(file_df)[1],
        exdir=dataset_dir)
    cel_files_gz <- list.files(dataset_dir, pattern=".*CEL.gz$",
        full.names=TRUE)
    for (f in cel_files_gz) GEOquery::gunzip(f, overwrite=TRUE)
    cel_files <- list.files(dataset_dir, pattern=".*CEL$",
        full.names=TRUE)
    rma_eset <- affy::justRMA(filenames=cel_files, cdfname=cdf)
    # drop control probes
    rma_eset <- rma_eset[!grepl("AFFX", rownames(rma_eset)), ]
    # format ENSG IDs
    rownames(rma_eset) <- gsub("_at$", "", rownames(rma_eset))
    # Fix sample names to match pData
    colnames(rma_eset) <- gsub("\\.CEL$", "", colnames(rma_eset))
    # add sample metadata
    phenoData(rma_eset) <- pData[colnames(rma_eset), ]
    # add feature metadata from gencode annotations
    feature_df <- data.table(gene_id=rownames(rma_eset))
    setkeyv(feature_df, "gene_id")
    fData <- gene_annots[feature_df, , on="gene_id"][!duplicated(gene_id), ]
    setDF(fData, rownames=fData$gene_id)
    featureData(rma_eset) <- as(fData, "AnnotatedDataFrame")
    # coerce to SummarizedExperiment
    se_type <- if (any(is.na(featureData(rma_eset)$start))) "SummarizedExperiment"
        else "RangedSummarizedExperiment"
    se_list[[ds]] <- as(rma_eset, se_type)
}

# Save datasets to disk
for (ds in names(se_list)) {
    qsave(se_list[[ds]],
        file=file.path(data_dir, paste0(ds, "_", class(se_list[[ds]])[1], ".qs")),
        nthread=getDTthreads()
    )
}