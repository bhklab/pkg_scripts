#!/bin/R

# -- Required
library(CoreGx)
library(BumpyMatrix)
library(SummarizedExperiment)

# -- Utilities
library(qs)
library(data.table)

## ==========================================
## ---- LongTable to gDR SummarizedExperiment
## ==========================================

# -- Load in the LongTable
LT <- qread('local_data/NCI_ALMANAC_LT.qs')

# ---- Convert assays to BumpyMatrix
# NOTE: We will need to pad the assays
# >if rows or columns are missing

#' Convert a LongTable assay into a BumpyMatrix object
#' 
#' @param LT `LongTable` with assay to convert into `BumpyMatrix`
#' @param assay `character(1)` A valid assay name in `LT`, as returned by 
#'     `assayNames(LT)`.
#' @param rows `character()` The rownames associated with the assay rowKey
#' @param cols `character()` The names associated with the assay colKey
#' @param sparse `logical(1)` Should the `BumpyMatrix` be sparse (i.e., is the
#'   assay sparse).
#' 
#' @return `BumpyMatrix` containing the data from `assay`.
#' 
#' @md
#' @importFrom CoreGx assay
#' @importFrom data.table data.table
.assayToBumpyMatrix <- function(LT, assay, rows, cols, sparse=TRUE) {
    assay_data <- assay(LT, assay, key=TRUE)
    assay_data[, rownames := rows[rowKey]]
    assay_data[, colnames := cols[colKey]]
    assay_data[, c('rowKey', 'colKey') := NULL]
    splitAsBumpyMatrix(assay_data[, -c('rownames', 'colnames')],
        row=assay_data$rownames, column=assay_data$colnames, sparse)
}

# ---- Convert LongTable to gDR SummarizedExperiment
#'
#' @param LT `LongTable` to convert to gDR `SummarizedExperiment` format.
#' @param assay_names `character()` Names to rename the assays to. These
#'   are assumed to be in the same order as `assayNames(LT)`.
#' 
#' @return `SummarizedExperiment` object with all assay from `LT` as 
#'   `BumpyMatrix`es.
#' 
#' @importFrom data.table setnames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom CoreGx rowData colData assayNames
.longTableToSummarizedExperiment <- function(LT, assay_names) {
    assay_list <- lapply(assayNames(LT), FUN=.assayToBumpyMatrix, 
        LT=LT, rows=rownames(LT), cols=colnames(LT))
    if (!missing(assay_names) && length(assay_names) == length(assayNames(LT))) 
        names(assay_list) <- assay_names
    SummarizedExperiment(
        assays=assay_list, rowData=rowData(LT), colData=colData(LT), 
            metadata=c(metadata(LT), list(.intern=as.list(getIntern(LT))))
    )
}

# -- build the SummarizedExperiment
SE <- .longTableToSummarizedExperiment(LT, 
    assay_names=c("Viability", "Metrics", "AssayMetadata"))