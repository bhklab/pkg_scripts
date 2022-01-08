library(qs)
library(CoreGx)
library(data.table)


# =================
# ---- Read in data


NCI <- qread("local_data/NCI_ALMANAC_LT.qs")

# -- Extract data
raw_data <- NCI$sensitivity


# =====================
# ---- Metadata storage


# -- rowData holds only treatment specification
rKeys <- c("drug1id", "drug2id")
rMeta <- cardinality(raw_data, rKeys)
rCols <- c(rKeys, rMeta)
rData <- unique(raw_data[, ..rCols])
rData[, rowKey := .I]

# -- colData holds only sample specification
cKeys <- c("cellid")
cMeta <- cardinality(raw_data, cKeys)
cCols <- c(cKeys, cMeta)
cData <- unique(raw_data[, ..cCols])
cData[, colKey := .I]


# ==================
# ---- Assay storage


# -- assays hold anything additional, such as doses, replicates, etc.
assay_raw <- copy(raw_data)
assay_raw[, rowKey := .GRP, by=.(drug1id, drug2id)]
assay_raw[, colKey := .GRP, by=.(cellid)]
assay_raw[, (c(rCols, cCols)) := NULL]


# ===================
# ---- Retrieve assay

# Join with metadata
a1 <- merge.data.table(assay_raw, rData, by="rowKey")
a1 <- merge.data.table(a1, cData, by="colKey")
setkeyv(a1, c("rowKey", "colKey"))


# =======================
# ---- Summarize an assay

# Drawback: lose ability to reconstruct the original data
# Benefit: subset with keys doesn't corrupt summaries, subset with assay does
a2 <- a1[, .SD, .SDcols=!c(rCols, cCols)][,
    c(lapply(.SD, mean), .(n=.N)),
    by=.(rowKey, colKey)
]


# =======================
# ---- Retrieve a summary

a3 <- merge.data.table(a2, rData, by="rowKey")
a3 <- merge.data.table(a3, cData, by="colKey")
setkeyv(a3, c("rowKey", "colKey"))


# ===============
# ---- Subsetting


# -- Subset the LongTable

# Get the matching keys
keep_rows <- rData[drug1id == "vemurafenib", ]$rowKey
keep_cols <- cData[cellid == "A498", ]$colKey

# Subset metadata
cData1 <- cData[colKey %in% keep_cols, ]
rData1 <- rData[rowKey %in% keep_rows, ]

# Subset assays
sub_expr <- quote(rowKey %in% keep_rows & colKey %in% keep_cols)
a1_1 <- a1[, .SD, .SDcols=!c(rCols, cCols)][i, , env=list(i=sub_expr)]
a2_1 <- a2[i, , env=list(i=sub_expr)]
a3_1 <- a3[, .SD, .SDcols=!c(rCols, cCols)][i, , env=list(i=sub_expr)]

# -- Subset by assay
## NOTE: This will corrupt summaries, do we want to allow that?
## Could drop the effected keys? That is probably a good idea

# Find the dimensions matching the criteria, in this case doses in a range
dose_range <- c(1.0e-10, 1.0e-8)
a1_2 <- a1[drug1dose %between% dose_range, ]
keep_dims <- lapply(a1_2[, .(rowKey, colKey)], FUN=unique)
i <- unique(a1_2[, .(rowKey, colKey)])

# Subset other assays
a2_2 <- a2[i, ]
a3_2 <- a3[i, ]

# Subset the metadata
cData2 <- cData[colKey %in% unique(keep_dims$colKey)]
rData2 <- rData[rowKey %in% unique(keep_dims$rowKey)]