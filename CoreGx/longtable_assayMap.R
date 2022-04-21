library(CoreGx)
library(data.table)
library(BiocParallel)

## ---- NCI-ALMANAC

# -- Configure script parameters
setDTthreads(2)
dataPath <- '../rawdata'
annotPath <- '../../Annotations'

# -- Read in data
exp_DT <- fread(file.path(dataPath,
    'NCI_ComboDrugGrowth_Nov2017.csv'))
drug_DT <- fread(file.path(dataPath,
    'NCI_compound_names.tsv'))
colnames(drug_DT) <- c('NSC', 'drug_name')

# -- Drop drugs with duplicated NSC ids
duplicatedDrugs <- drug_DT[duplicated(NSC) | duplicated(drug_name), ]
drug_DT <- drug_DT[!(duplicated(NSC) | duplicated(drug_name)), ]

# -- Join to map to names of drug1
#    and drug2 in combo
nci_DT <- merge.data.table(exp_DT, drug_DT, by.x='NSC1', by.y='NSC')
setnames(nci_DT,
    old=c('drug_name', 'CONC1'),
    new=c('drug1id', 'drug1dose'),
    skip_absent=TRUE)

nci_DT <- merge.data.table(nci_DT, drug_DT, by.x='NSC2', by.y='NSC')
setnames(nci_DT,
    old=c('drug_name', 'CONC2', 'CELLNAME'),
    new=c('drug2id', 'drug2dose', "cellid"),
    skip_absent=TRUE)


# -- Define a data mapper class
dataMapperLT <- TREDataMapper(rawdata=nci_DT)

# -- Guess the mappings to different identifiers
groups <- list(
    colDataMap=c("cellid"),
    rowDataMap=c("drug1id", "drug2id", "drug1dose", "drug2dose"),
    assayMap=c("cellid", "drug1id", "drug2id", "drug1dose", "drug2dose")
)

## TODO:: Implement the faster version of this using a matrix
guess <- guessMapping(dataMapperLT, groups=groups, subset=TRUE)

assayMap_ <- list(
    sensitivity=list(
        id_columns=groups$assayMap,
        mapped_columns=c("CONCINDEX2", "SAMPLE1", "SAMPLE2",
            viability='PERCENTGROWTH', 'PERCENTGROWTHNOTZ',
            'EXPECTEDGROWTH', "TESTVALUE", "CONTROLVALUE", "TZVALUE")
    ),
    profiles=list(
        id_columns=groups$assayMap,
        mapped_columns=c('SCORE')
    ),
    assay_metadata=list(
        id_columns=groups$assayMap,
        mapped_columns=c('COMBODRUGSEQ', 'STUDY', 'TESTDATE', 'PLATE')
    )
)

## -- Define mappings for row and column annotations
rData <- unique(nci_DT[, .SD, .SDcols=unlist(guess$rowDataMap)])
rData[, rowKey := .I, by=c(guess$rowDataMap$id_columns)]
setkeyv(rData, "rowKey")

cData <- unique(nci_DT[, .SD, .SDcols=unlist(guess$colDataMap)])
cData[, colKey := .I, by=c(guess$colDataMap$id_columns)]
setkeyv(cData, "colKey")

## -- Extract metadata
metadata_ <- unique(nci_DT[, .SD, .SDcols=guess$metadata[[2]]])

## -- Build assays
aKeys <- groups$assayMap

nci_DT[, rowKey := .GRP, by=c(guess$rowDataMap$id_columns)]
nci_DT[, colKey := .GRP, by=c(guess$colDataMap$id_columns)]
nci_DT[, ogKey := .GRP, keyby=.(rowKey, colKey)]

a1 <- unique(nci_DT[,
    .SD, .SDcols=c("ogKey", assayMap_$sensitivity$mapped_columns)])
a2 <- unique(nci_DT[,
    .SD, .SDcols=c("ogKey", assayMap_$profiles$mapped_columns)])
a3 <- unique(nci_DT[,
    .SD, .SDcols=c("ogKey", assayMap_$assay_metadata$mapped_columns)])
assays <- list(a1, a2, a3)
for (a in assays) setkeyv(a, "ogKey")

## -- Coerce back to data.table of raw data
aMap <- unique(nci_DT[, .(rowKey, colKey, ogKey)])
setkeyv(aMap, c("rowKey", "colKey", "ogKey"))

obj_ <- c(list(rData, cData), assays)
# join everything with aMap
obj_ <- lapply(obj_, merge, y=aMap)
# set key to all keys
for (o_ in obj_) setkeyv(o_, c("ogKey", "rowKey", "colKey"))
# reduce merge
dt_ <- Reduce(f=merge, obj_)
dt_ <- cbind(dt_, as.data.table(metadata_))


## -- metaConstruct works
## FIXME:: Why is this assignment slow?
rowDataMap(dataMapperLT) <- guess$rowDataMap
colDataMap(dataMapperLT) <- guess$colDataMap
assayMap(dataMapperLT) <- assayMap_
metadataMap(dataMapperLT) <- list(
    experiment_metadata=guess$metadata$mapped_columns
)
