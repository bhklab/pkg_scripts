library(data.table)
library(SummarizedExperiment)
library(downloader)
library(rvest)
library(qs)
library(readxl)
library(BiocParallel)
library(progress)

# get annotation package for Human Exon 1.0 ST transcript clusters
# NOTE: probset_id in the TARGET OS annotations is actually transcript cluster id!
if (!require("huex10sttranscriptcluster.db"))
    BiocManager::install("huex10sttranscriptcluster.db")
library(huex10sttranscriptcluster.db)

# load array annotation datasets
huex10sttranscriptcluster()
probeset_to_ensembl <- as.data.table(huex10sttranscriptclusterENSEMBL)
probeset_to_entrez <- as.data.table(huex10sttranscriptclusterENTREZID)

# paste together ids for probesets mapping to multiple genes
## NOTE: entrez tables have no multimaps
probeset_to_ensembl_no_dups <- probeset_to_ensembl[,
    .(ensembl_id=paste0(ensembl_id, collapse="|")),
    by=probe_id
]

# -- script configuration parameters
data_dir <- "local_data/TARGET_OS/gene_expression_array"
metadata_dir <- "metadata"

url <- "https://target-data.nci.nih.gov/Public/OS/gene_expression_array"

# -- helper function definitions
grep_directory_names <- function(x) {
    grep(".*/$", x, value=TRUE)
}

grep_file_names <- function(x, extensions="[[:alnum:]]{2,5}$") {
    checkmate::assert_character(extensions, max.len=1, min.chars=1)
    grep(paste0(".*\\.", extensions), x, value=TRUE)
}

get_remote_table <- function(url) {
    read_html(url) |>
    html_elements("body") |>
    html_table()
}

#' Recursively find file URLs from URL with an embedded HTML table (a remote directory)
#'
#' @param url `character(1)` A valid URL to scrape file data from. It is assumed
#'   that the returned HTML has a table in it which indicates the remote
#'   directory contents.
#' @param column `character(1)` Name of the column in the returned HTML table
#'   to match files and directories from. The `url` is automatically prepended
#'   to the values in the column, so they should be relative paths to other files
#'   in the remote directory.
#' @param extensions `character()` vector with one of more file extensions to
#'   to scrape from `url`. This should be the file extension only, with no dot.
#'   It could also be a valid regex expression. Please note that all values will
#'   be appended with "$" to match on the end of files and collapsed together
#'   with the "|" regex operator. The default matches any alphanumeric file
#'   extensions between two and five character long.
#'
#' @return `character()` vector of remote file URLs to download from.
#'
#' @importFrom checkmate assert_character
#' @export
find_remote_files_recursive <- function(url, column="Name", extensions="[[:alnum:]]{2,5}") {
    # input validation
    checkmate::assert_character(url, min.chars=5, max.len=1)
    checkmate::assert_character(column, min.chars=1, max.len=1)
    checkmate::assert_character(extensions, min.chars=2, min.len=1)
    exts <- paste0(paste0(gsub("\\.", "", extensions), "$"), collapse="|")
    # extract all HTML tables
    tables <- get_remote_table(url)
    # extract directory and file names
    directories <- gsub("/$", "",
        unlist(lapply(tables, FUN=\(x) grep_directory_names(x[[column]]))))
    files <- unlist(lapply(tables, FUN=\(x) grep_file_names(x[[column]], exts)))
    # append files and directories to current url
    directories <- file.path(url, directories)
    files <- file.path(url, files)
    if (length(directories > 0)) {
        # directories become new URL for recursive calls
        new_files <- unlist(lapply(directories, FUN=find_remote_files_recursive))
        files <- c(files, new_files)
    } else {
        return(files)
    }
}

# -- download the raw data
remote_files <- find_remote_files_recursive(url)
file_paths <- gsub(url, data_dir, remote_files)

for (i in seq_along(remote_files)) {
    if (!dir.exists(dirname(file_paths[i])))
        dir.create(dirname(file_paths[i]), recursive=TRUE)
    if (!file.exists(file_paths[i]))
        downloader::download(remote_files[i], destfile=file_paths[i])
}

# -- Read in the sample metadata
list.files(data_dir)
metadata_path <- list.files(file.path(data_dir, "METADATA"), full.names=TRUE)[3]
col_data <- fread(metadata_path)

# -- Read in gene expression array
gene_expr_files <- list.files(
    file.path(data_dir, "L3"),
    pattern="gene",
    full.names=TRUE
)
names(gene_expr_files) <- gene_expr_files
gene_expr <- lapply(gene_expr_files, FUN=fread)[[1]]

# -- Map the probset_id column to associated ensembl genes

# -- Parse the gene names into valid regex queries
gene_expr[,
    gene_id := unlist(lapply(gene_assignment_final,
        FUN=\(x) paste0(unique(strsplit(x, " // ")[[1]]), collapse="|")))
]

expr_mat <- as.matrix(gene_expr[, -c("gene_id", "gene_assignment_final")])
rownames(expr_mat) <- gene_expr$probeset_id

# -- Subset the sample metadata to only included samples
setkeyv(col_data, "Array Data File")
subColData <- col_data[colnames(expr_mat), ]
subColData <- subColData[, first(.SD), by=`Array Data File`]

# -- Make SummarizedExperiment
target_os <- SummarizedExperiment(
    assays=list(rma=expr_mat),
    colData=subColData,
    rowData=gene_expr[, .SD, .SDcols=!patterns(".CEL")]
)

# -- Read and merge clinical metadata
clinical_path <- file.path("local_data", "TARGET_OS", "clinical")
clinical_files <- list.files(clinical_path, pattern="xlsx", full.names=TRUE)
clinical_meta <- lapply(clinical_files, FUN=read_excel)

metadata(target_os)$clinical_data_elements <- as.data.frame(clinical_meta[[1]])

# match the colnames to TARGET USI
colnames(target_os) <- colData(target_os)$`Source Name`
colData(target_os)$`TARGET USI` <- colnames(target_os)

clinical_dates <- gsub(".*_|.xlsx", "", clinical_files[-1])
clinical_dts <- lapply(clinical_meta[-1], as.data.table)
clinical_dts <- lapply(clinical_dts,
    FUN=\(x, cols) x[`TARGET USI` %in% cols, ],
    cols=colnames(target_os)
)
for (i in seq_along(clinical_dts)) {
    colnames(clinical_dts[[i]])[-1] <- paste0(clinical_dates[i], ".",
        colnames(clinical_dts[[i]])[-1])
    setkeyv(clinical_dts[[i]], "TARGET USI")
}
merge_all <- function(x, y) merge(x, y, all=TRUE, by="TARGET USI")
clinical_meta <- Reduce(merge_all, clinical_dts[-2])

colData <- as.data.table(as.list(colData(target_os)))
colData <- merge(colData, clinical_meta, by="TARGET USI", all.x=TRUE)
colData(target_os) <- DataFrame(colData, check.names=FALSE)
colnames(target_os) <- colData(target_os)$`TARGET USI`


# -- Add gene annotations to rowData

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
gene_annots[, gene_id_no_ver := gsub("\\..*$", "", gene_id)]
setnames(probeset_to_ensembl_no_dups,  # clean up names to indicate source of annotation
    c("probe_id", "ensembl_id"),
    c("probeset_id", "huex10sttranscriptcluster.ensembl_id")
)
probeset_to_gencode <- merge.data.table(
    probeset_to_ensembl_no_dups,
    gene_annots,
    by.x="huex10sttranscriptcluster.ensembl_id",
    by.y="gene_id_no_ver",
    all.x=TRUE
)
# Add the multimapped genes as well
multi_gene_patterns <- probeset_to_ensembl_no_dups[
    huex10sttranscriptcluster.ensembl_id %like% "\\|",
][["huex10sttranscriptcluster.ensembl_id"]]
multi_gene_matches <- vector("list", length(multi_gene_patterns)) |>
    setNames(multi_gene_patterns)
pb <- progress_bar$new(total=length(multi_gene_matches))
for (pattern in multi_gene_patterns) {
    pb$tick()
    multi_gene_matches[[pattern]] <- gene_annots[
        gene_id_no_ver %like% pattern,
    ]
}
multi_gene_dt <- rbindlist(multi_gene_matches, idcol=TRUE)
# select only the first gene where multimapping occurs
multi_gene_no_multi <- multi_gene_dt[, first(.SD), by=.id]
multi_gene_mapped <- merge.data.table(
    probeset_to_ensembl_no_dups,
    multi_gene_no_multi,
    by.x="huex10sttranscriptcluster.ensembl_id",
    by.y=".id"
)
# Merge clean mapped with multimapped gencode annotations
probeset_to_gencode <- rbind(
    probeset_to_gencode,
    multi_gene_mapped[, colnames(probeset_to_gencode), with=FALSE]
)

# Attach annotations to rowData of the TARGET OS SummarizedExperiment
rdata <- as.data.table(rowData(target_os))[, gene_id := NULL]
rdata[, probeset_id := as.character(probeset_id)]
rdata <- merge.data.table(rdata, probeset_to_gencode, by="probeset_id",
    all.x=TRUE, sort=FALSE)
setnames(rdata, "gene_assignment_final", "target_transcript_id")
# Drop duplicated probeset ids by selecting those with valid gene_ids or
#   if there are none, select the first occurence
rdata <- rdata[,
    if (.N > 1)
        .SD[, if (any(!is.na(gene_name))) .SD[!is.na(gene_name), ][1, ] else first(.SD)]
    else .SD,
    by=probeset_id
]

# Sanity check that the gene symbols match the original TARGET annotations
mismatches <- list()
num_valid_genes <- nrow(rdata[!is.na(gene_name)])
pb <- progress_bar$new(total=num_valid_genes)
for (i in seq_len(num_valid_genes)) {
    pb$tick()
    row_res <- rdata[!is.na(gene_name)][
        i,
        if (target_transcript_id %like% gene_name) TRUE else .SD
    ]
    if (!isTRUE(row_res)) {
        mismatches <- append(mismatches, list(row_res))
    }
}

# Assign the feature annotations back to the object
setkeyv(rdata, "probeset_id")
rowData(target_os) <- rdata[rownames(target_os), ]

# -- Save file to disk
qsave(target_os, file=file.path(data_dir, "target_os_micro_se.qs"), nthread=8)