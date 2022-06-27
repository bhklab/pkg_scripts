library(data.table)
library(rtracklayer)

# no timeouts on downloads
options(timeout=Inf)

# -- Download the data files

datadir <- ".local_data"
gencodev40_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"
gencodev40_uniprot_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.metadata.SwissProt.gz"


urls <- c(gencodev40_url, gencodev40_uniprot_url)
file_paths <- vector("character", length(urls))
for (i in seq_along(urls)) {
    file_name <- basename(urls[i])
    file_paths[i] <- file.path(datadir, file_name)
    download.file(
        urls[i],
        destfile=file_paths[i]
    )
}
gencodev40 <- import(gzfile(file_paths[1]))
gencode_dt <- as.data.table(gencodev40)

gencode_files <- lapply(file_paths[-1], fread, header=FALSE, sep="\t")
uniprot_annot <- gencode_files[[1]]
colnames(uniprot_annot) <- c("transcript_id", "uniprot_id", "uniprot_id_ver")

gencode_uniprot <- merge.data.table(
    gencode_dt,
    uniprot_annot,
    by="transcript_id"
)

gencodev40_valid_uniprot <- unique(gencode_uniprot[!is.na(uniprot_id)])
fwrite(
    gencodev40_valid_uniprot, 
    file=file.path(datadir, "gencodev40_to_uniprot.csv")
)