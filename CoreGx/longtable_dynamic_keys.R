library(qs)
library(CoreGx)
library(data.table)


# =================
# ---- Read in data


NCI <- qread("local_data/NCI_ALMANAC_LT.qs")

# -- Extract components
rData <- rowData(NCI)
cData <- colData(NCI)
assay1 <- assay(NCI, "sensitivity")


# ==========
# ---- Joins


# Keys will be attributes attached to assays of the object
rkey1 <- c("drug1id", "drug2id", "drug1dose", "drug2dose")
ckey1 <- c("cellid")
rkey2 <- c("drug1id", "drug2id")  # Summarize over dose

# -- Retrieving an assay with metadata
# copy the metadata
rd1 <- copy(rData)
cd1 <- copy(cData)

# compute a key on the fly
rd1[, rowKey := .GRP, by=rkey1]
cd1[, colKey := .GRP, by=ckey1]

a1_joined <- merge.data.table(assay1, rd1, by="rowKey")
a1_joined <- merge.data.table(a1_joined, cd1, by="colKey")
a1_joined[, (c("rowKey", "colKey")) := NULL]

# -- Retrieving an assay with metadata for a summary
rd2 <- copy(rd1)
rd2[, rowKey := .GRP, by=rkey2]

# Make a summary assay (how to encode this as a function?)
assay2 <- a1_joined[, .(mean_viability=mean(viability)), by=c(rkey2, ckey1)]
assay2[, rowKey := .GRP, by=rkey2]
assay2[, colKey := .GRP, by=ckey1]
assay2[, (c(rkey2, ckey1)) := NULL]

# Reattach metadata (requires cartesian join!)
a2_joined <- merge.data.table(rd2, assay2, by="rowKey", allow.cartesian=TRUE)
a2_joined <- merge.data.table(cd1, assay2, by="colKey", allow.cartesian=TRUE)

# Why are there more rows?
# - Losing information about which doses map to which cells, effectively
#   dropping a join table needed to reconstruct the unsummarized data

# ============
# ---- Subsets





# ================
# ---- Aggregation



# ==============
# ---- Modelling