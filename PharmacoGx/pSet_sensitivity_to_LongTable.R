#!/bin/R

# -- Required
library(PharmacoGx)
library(CoreGx)
library(gDR)
library(data.table)
library(qs)
library(BumpyMatrix)


# -- Download some data
dataset <- "gCSI"
data_dir <- "local_data"
availPSets <- as.data.table(availablePSets())
pset_name <- availPSets[`Dataset Name` == dataset, ]$`PSet Name`

if (!file.exists(file.path(data_dir, dataset))) {
    assign(dataset,
        downloadPSet(
            pset_name,
            saveDir=data_dir
        )
    )
} else {
    assign(dataset, readRDS(file.path(data_dir, pset_name)))
}

## -- Convert the sensitivity slot to a TreatmentResponseExperiment
tre <- CoreGx:::.sensitivityToLongTable(get(dataset))

## -- Convert the TreatmentResponseExperiment to a long format data.table
## TODO:: make as.data.table work
dt <- as(tre, "data.table")
## FIXME:: why is TrtDuration is getting dropped inside .sensitivityToLongTable?
dt$duration <- 72
df_ <- copy(dt[
        cellid %in% unique(cellid)[1:5] & drug1id %in% unique(drug1id)[1:5],
        .(cellid, doublingtime, drug1id, drug1dose, viability, duration)
])

# -- Free up memory
rm(list=c("tre", dataset))
gc()

## -- Set environment identifiers to map from our columns to gDR columns
gDRutils::reset_env_identifiers()
default_ids <- gDRutils::get_env_identifiers()
## TODO:: helper function to guess identifier mappings
pgx_to_gdr_ids <- list(
    cellline="cellid",
    cellline_name="cellid",
    cellline_ref_div_time="doublingtime",
    drug="drug1id",
    drugname="drug1id",
    concentration="drug1dose",
    duration="duration"
)
unmodified_ids <- default_ids[setdiff(names(default_ids), names(pgx_to_gdr_ids))]
pgx_to_gdr_ids <- c(pgx_to_gdr_ids, unmodified_ids)
## FIXME:: Setting identifiers drops defaults?
for (i in seq_along(pgx_to_gdr_ids))
    gDRutils::set_env_identifier(
        k=names(pgx_to_gdr_ids)[i],
        v=pgx_to_gdr_ids[[i]]
    )

## -- Create the gDR SummarizedExperiment
se <- create_SE(
    df_,
    readout="viability",
    nested_keys=get_env_identifiers("concentration")
)


# ====================
# ==== DEPRECATED ====
# ====================


# ---- Convert assays to BumpyMatrix
# NOTE: We will need to pad the assays
# >if rows or columns are missing

#' Convert a LongTable assay into a BumpyMatrix object
#'
#' @param TRE `TreatmentResponseExperiment` with assay to convert into
#'     `BumpyMatrix`
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
.assayToBumpyMatrix <- function(TRE, assay, rows, cols, sparse=TRUE) {
    assay_data <- assay(TRE, assay, key=TRUE)
    assay_data[, rownames := rows[rowKey]]
    assay_data[, colnames := cols[colKey]]
    assay_data[, c('rowKey', 'colKey') := NULL]
    splitAsBumpyMatrix(assay_data[, -c('rownames', 'colnames')],
        row=assay_data$rownames, column=assay_data$colnames, sparse)
}

#' Convert LongTable to gDR Style SummarizedExperiment
#'
#' @param TRE `LongTable` to convert to gDR `SummarizedExperiment` format.
#' @param assay_names `character()` Names to rename the assays to. These
#'   are assumed to be in the same order as `assayNames(TRE)`.
#' 
#' @return `SummarizedExperiment` object with all assay from `TRE` as 
#'   `BumpyMatrix`es.
#' 
#' @importFrom data.table setnames
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom CoreGx rowData colData assayNames
.trExperimentToSummarizedExperiment <- function(TRE, assay_names) {
    assay_list <- lapply(assayNames(TRE), FUN=.assayToBumpyMatrix, 
        TRE=TRE, rows='drug1id', cols='cellid')
    if (!missing(assay_names) && length(assay_names) == length(assayNames(TRE))) 
        names(assay_list) <- assay_names
    SummarizedExperiment(
        assays=assay_list, rowData=rowData(TRE), colData=colData(TRE), 
            metadata=c(metadata(TRE), list(.intern=as.list(getIntern(TRE))))
    )
}


se <- .longTableToSummarizedExperiment(tre,
    assay_names=c('Viability', 'Metrics', 'AssayMetadata'))

qsave(se, file=file.path(datadir, paste0(dataset, "_gDR.qs")))