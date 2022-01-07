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
a1[, (c("rowKey", "colKey")) := NULL]


# =======================
# ---- Summarize an assay

# Drawback: lose ability to reconstruct the original data
# Benefit: subset with keys doesn't corrupt summaries, subset with assay does
a2 <- a1[, !c("PANEL", "SCREENER")][,
    c(lapply(.SD, mean), .(n=.N)),
    by=c(rKeys, cKeys)
]


# =======================
# ---- Retrieve a summary

a3 <- merge.data.table(a2, rData, by=rKeys)
a3 <- merge.data.table(a3, cData, by=cKeys)
a3[, (c("rowKey", "colKey")) := NULL]


# ===============
# ---- Subsetting


# -- Subset the LongTable

# Get the matching keys
keep_rows <- rData[drug1id == "vemurafenib", ]$rowKey
keep_cols <- cData[cellid == "A498", ]$colKey

# Subset metadata
cData1 <- cData[colKey %in% keep_cols, ]
rData1 <- rData[rowKey %in% keep_rows]

# Subset assays
sub_expr <- quote(rowKey %in% keep_rows & colKey %in% keep_cols)
a11 <- a1[i, , env=list(i=sub_expr)]
a22 <- a2[rowKey %in% keep_rows & colKey %in% keep_cols, ]
a33 <- a3[]

# -- Subset by assay

## NOTE: This will corrupt summaries, do we want to allow that?
