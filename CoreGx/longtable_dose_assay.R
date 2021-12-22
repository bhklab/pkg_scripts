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
rData <- unique(raw_data[, ..rKeys])
rData[, rowKey := .I]

# -- colData holds only sample specification
cKeys <- c("cellid")
cData <- unique(raw_data[, ..cKeys])
cData[, colKey := .I]

# ==================
# ---- Assay storage


# -- assays hold anything additional, such as doses, replicates, etc.
assay_raw <- copy(raw_data)
assay_raw[, rowKey := .GRP, by=.(drug1id, drug2id)]
assay_raw[, colKey := .GRP, by=.(cellid)]


# ===================
# ---- Retrieve assay

# Join with metadata
a1 <- merge.data.table(assay_raw, rData, by="rowKey")
a1 <- merge.data.table(a1, cData, by="colKey")


# =======================
# ---- Summarize an assay

# Drawback: lose ability to reconstruct the original data
a2 <- a1[, c(lapply(.SD, mean), .(n=.N)), by=c(rKeys, cKeys)]
a2[, (c(rKeys, cKeys)) := NULL]