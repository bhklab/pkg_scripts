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

# ---- Capture metadata for SummarizedExperiment
rowDat <- rowData(LT)
rowDat$rownames <- rownames(LT)
colDat <- colData(LT)
colDat$rownames <- colnames(LT)

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
#' @param ... Fallthrough arguments to splitAsBumpyMatrix.
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

# -- Convert assays to `BumpyMatrix`es
assay_list <- lapply(assayNames(LT), .assayToBumpyMatrix, LT=LT, 
    rows=rowDat$rownames, cols=colDat$rownames)

# -- build the SummarizedExperiment
SE <- SummarizedExperiment(
    assays=list(Viability=viab_BM, Profiles=prof_BM, 
        AssayMetadata=meta_BM),
    rowData=rowDat, colData=colDat, metadata=metadata(LT))