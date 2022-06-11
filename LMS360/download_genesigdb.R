library(rvest)
library(jsonlite)
library(data.table)
library(future)
library(future.apply)
library(R.utils)

# Script parameters
data_dir <- "local_data"

# Get the GeneSigDB metadata
harmonizome_url <- "https://maayanlab.cloud/Harmonizome"
genesigdb_json_url <- paste0(harmonizome_url, "/api/1.0/dataset/GeneSigDB+Published+Gene+Signatures")
genesigdb_list <- fromJSON(genesigdb_json_url)

# Extract the table of gene signature API urls
geneset_df <- as.data.table(genesigdb_list$geneSets)
geneset_df[, geneset_url := paste0(harmonizome_url, href)]

# Set our execution plan to use forked parallization (doesn't work on Windows)
future::plan("multicore")

# Do asynchronous calls to fetch the signatures from the API
genesets <- future_lapply(geneset_df$gene, FUN=fromJSON)

# Parse the signature genesets into a long table'
geneset_dfs <- future_lapply(genesets, FUN=\(x) {
    gene_df <- x$associations$gene
    gene_df$SigID <- x$attribute$name
    gene_df
})
signature_df <- rbindlist(geneset_dfs)
signature_df[, href := NULL]

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

# Attach gene annotations
signature_df <- merge.data.table(signature_df,
    gene_annots, by.x="symbol", by.y="gene_name",
    all.x=TRUE
)

# Get index with signature descriptions from GeneSigDB GitHub repo
signature_metadata_url <- "https://github.com/aedin/GeneSigDB/raw/master/data/GeneSigDB-Table%201.csv"
signature_metadata <- fread(signature_metadata_url)

# Drop useless columns
signature_metadata <- signature_metadata[, .SD, .SDcols=!patterns("Column\\d+|V\\d+")]

signature_df <- merge.data.table(signature_df, signature_metadata, by="SigID",
    all.x=TRUE)

fwrite(signature_df, file=.file.path(data_dir, "genesigdb_all_signature.csv"))

sarcoma_sigs <- signature_df[Tissue == "Sarcoma", ]