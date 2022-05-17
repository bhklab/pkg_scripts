library(data.table)
library(SummarizedExperiment)
library(downloader)
library(rvest)
library(qs)
library(readxl)

data_dir <- "local_data/TARGET_OS/rnaseq"
metadata_dir <- "metadata"

url <- "https://target-data.nci.nih.gov/Public/OS/mRNA-seq"

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
    downloader::download(remote_files[i], destfile=file_paths[i])
}


# -- Read in the sample metadata
list.files(data_dir)
metadata_path <- list.files(file.path(data_dir, "METADATA"), full.names=TRUE)[3]
col_data <- fread(metadata_path)

# -- Read in gene expression
gene_expr_files <- list.files(
    file.path(data_dir, "L3", "expression", "NCI-Meltzer"),
    pattern="gene",
    full.names=TRUE
)
names(gene_expr_files) <- gene_expr_files
gene_expr <- lapply(gene_expr_files, FUN=fread)
gene_expr_long <- rbindlist(gene_expr, idcol="file")

# -- Parse the long table into a gene by sample matrix
gene_expr_long[, file := basename(file)]
setnames(gene_expr_long, "Name", "gene_id")
gene_expr_wide <- dcast(gene_expr_long, gene_id ~ file, value.var="TPM")

expr_mat <- as.matrix(gene_expr_wide[, -c("gene_id")])
rownames(expr_mat) <- gene_expr_wide$gene_id

# -- Subset the sampe metadata to only included samples
setkeyv(col_data, "Derived Array Data File")
subColData <- col_data[colnames(expr_mat), ]
subColData <- subColData[, first(.SD), by=`Derived Array Data File`]

# -- Make SummarizedExperiment
target_os <- SummarizedExperiment(assays=list(tpm=expr_mat), colData=subColData)

# -- Download clinical metadata
url <- "https://target-data.nci.nih.gov/Public/OS/clinical/harmonized"
remote_files <- find_remote_files_recursive(url)
clinical_path <- file.path(data_dir, "clinical")
local_files <- gsub(url, clincal_path, remote_files)

for (i in seq_along(remote_files)) {
    if (!dir.exitst(dirname(local_files[i])))
        dir.create(dirname(local_files[i]), recursive=TRUE)
    downloader::download(remote_files[i], destfile=local_files[i])
}

# -- Read and merge clinical metadata
clinical_files <- list.files(clinical_path, pattern="xlsx", full.names=TRUE)
clinical_meta <- lapply(clinical_files, FUN=read_excel)

metadata(target_os)$clinical_data_elements <- as.data.frame(clinical_meta[[1]])

# match the colnames to TARGET USI
colnames(target_os) <- gsub("-[^-]*-[^-]*\\.gene.quantification.txt", "",
    colnames(target_os))
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

qsave(target_os, file=file.path(data_dir, "target_os_se.qs"), nthread=8)