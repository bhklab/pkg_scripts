library(data.table)
library(qs)
library(PharmacoGx)

# script parameters
nthread <- 10

# read in data
sarcsets <- qread("local_data/sarcsets.qs", nthread=nthread)  # From download_pgx_sarcoma_data.R
cell_annots <- fread("local_data/cell_annotation_all.csv")  # From BHKLAB-Pachyderm/Annotations
dimitrios_annots <- fread("metadata/NCI_Sarcoma_cellInfo_bonevssoft_sarcoma_drdimitrios.csv")

# add tissueid to NCI_Sarcoma cellInfo
sarc_cInfo <- cellInfo(sarcsets$NCI_Sarcoma)
setDT(sarc_cInfo, keep.rownames=TRUE)

cellid_to_tissueid <- cell_annots[,
    .(cellid=unique.cellid, tissueid=unique.tissueid)
]

cInfo <- merge.data.table(sarc_cInfo, cellid_to_tissueid,
    by="cellid", all.x=TRUE)
setcolorder(cInfo, c("rn", "cellid", "tissueid"))
cellInfo(sarcsets$NCI_Sarcoma) <- as.data.frame(cInfo[, -"rn"],
    row.names=cInfo[["rn"]])

cellInfo(sarcsets$NCI_Sarcoma) <- merge(cellInfo(sarcsets$NCI_Sarcoma),
    dimitrios_annots[, .(cellid, dimitrios.soft_vs_bone=`Soft vs Bone`)],
    by="cellid", all.x=TRUE)
# fix rownames
rownames(cellInfo(sarcsets$NCI_Sarcoma)) <- cellInfo(sarcsets$NCI_Sarcoma)$cellid

# write back to disk
qsave(sarcsets, file="local_data/sarcsets.qs", nthread=nthread)