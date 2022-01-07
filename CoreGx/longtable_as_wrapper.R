library(qs)
library(CoreGx)

# ---- Load data
NCI <- qread("local_data/NCI_ALMANAC_LT.qs")
raw <- as(NCI, "data.table")


# ======================
# ---- Set wrapper class

setClass("DataTableExperiment", slots=list(
    .data="data.table",
    .assays="list",
    .row="list",
    .cols="list",
    metadata="list"
))

dt_array <- new("DataTableArray",
    .data=raw,

)

setMethod("assay", signature())