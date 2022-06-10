library(CoreGx)
library(BiocParallel)
library(data.table)

## ---- NCI-ALMANAC

# -- Configure script parameters
setDTthreads(2)
dataPath <- '.local_data'

# -- Read in data
exp_DT <- fread(file.path(dataPath,
    'NCI_ComboDrugGrowth_Nov2017.csv'))
drug_DT <- fread(file.path(dataPath,
    'NCI_compound_names.tsv'))
colnames(drug_DT) <- c('NSC', 'drug_name')

# -- Check drugs with duplicated NSC ids
duplicatedDrugs <- drug_DT[duplicated(NSC) | duplicated(drug_name), ]
# drug_DT <- drug_DT[!(duplicated(NSC) | duplicated(drug_name)), ]
## shouldn't remove duplicated drug names (i.e. the same drug with multiple NSC mapped to)
## This will cause 753082 vemurafenib to be removed from drug_DT but it is in exp_DT

# -- Combine NSC ids with multiple drug names into a single value
drug_DT <- drug_DT[, .(drug_name=paste0(drug_name, collapse="|")), by=NSC]

# -- Handle drugs with duplicated NSC ids
drug_DT[duplicated(drug_name), drug_name := paste0(drug_name, "2")]

# -- Join to map to names of drug1
#    and drug2 in combo
nci_DT <- merge.data.table(exp_DT, drug_DT, by.x='NSC1', by.y='NSC')
if (dim(exp_DT)[1] != dim(nci_DT)[1]) {
    stop("Rows lost in exp_DT")
}
setnames(nci_DT,
    old=c('drug_name', 'CONC1'),
    new=c('drug1id', 'drug1dose'),
    skip_absent=TRUE)

# nci_DT <- merge.data.table(nci_DT, drug_DT, by.x='NSC2', by.y='NSC')
## This line results in loss of monothreapy data,
## leading to only 2810343 observations
## We use left join instead
nci_DT <- merge.data.table(nci_DT, drug_DT, by.x='NSC2', by.y='NSC', all.x = TRUE)
if (dim(exp_DT)[1] != dim(nci_DT)[1]) {
    exp_DT[!(NSC2 %in% nci_DT[, NSC2]),]
    stop("Rows lost in exp_DT")
}
setnames(nci_DT,
    old=c('drug_name', 'CONC2', 'CELLNAME'),
    new=c('drug2id', 'drug2dose', "cellid"),
    skip_absent=TRUE)

# -- Define a data mapper class
dataMapperTRE <- TREDataMapper(rawdata=nci_DT)

# -- Guess the mappings to different identifiers
groups <- list(
    colDataMap=c("cellid"),
    rowDataMap=c("drug1id", "drug2id", "drug1dose", "drug2dose"),
    assayMap=c("cellid", "drug1id", "drug2id", "drug1dose", "drug2dose")
)

## TODO:: Implement the faster version of this using a matrix
#guess <- guessMapping(dataMapperTRE, groups=groups, subset=TRUE)
# bench::system_time({
#     guess <- guessMapping(dataMapperTRE, groups=groups, subset=TRUE)
# })
#process    real
#  4.45m   4.36m
## This is significantly slower than the case without monotherapy data (13s)
## Oddly, using only 2 threads took less time than using more threads:
#process    real
#  3.71m   3.68m
## Chris: yea with data.table using more threads isn't always faster, it is
## pretty hardware dependent. In my experience, anything over 8 threads ends up
## having too much parallelization overhead to see any speed ups

# -- handle case where assay keys are insufficient to uniquely identify rows
# if (length(guess$unmapped) > 0) {
    raw <- rawdata(dataMapperTRE)
    raw[, replicate_id := seq_len(.N), by=c(groups$assayMap)]
    nci_DT[, replicate_id := seq_len(.N), by=c(groups$assayMap)]
    groups$rowDataMap <- c(groups$rowDataMap, "replicate_id")
    groups$assayMap <- c(groups$assayMap, "replicate_id")
    rawdata(dataMapperTRE) <- raw
    bench::system_time({
        guess <- guessMapping(dataMapperTRE, groups=groups, subset=TRUE)
    })
# }


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
dt_ <- Reduce(f=merge, obj_) ## R barks on this line, results in the following error message:

dt_ <- cbind(dt_, as.data.table(metadata_))

rm(list=setdiff(ls(), c("dataMapperTRE", "guess", "assayMap_"))); gc()

## -- metaConstruct works
## FIXME:: Why is this assignment slow?
bench::system_time({ rowDataMap(dataMapperTRE) <- guess$rowDataMap })
bench::system_time({ colDataMap(dataMapperTRE) <- guess$colDataMap })
bench::system_time({ assayMap(dataMapperTRE) <- assayMap_ })
bench::system_time({
    metadataMap(dataMapperTRE) <- list(
        experiment_metadata=guess$metadata$mapped_columns
    )
})

bench::system_time({ tre <- metaConstruct(dataMapperTRE) })

## -- LongTable demo

## rowData (treatments)

# view the rowData
rowData(tre)

# see the row identifier column names
rowIDs(tre)
# see the row metadata column names
rowMeta(tre)

## colData (samples)

# view the coldata
colData(tre)

# see the column identifier column names
colIDs(tre)
# see the column metdata colum names
colMeta(tre)

## viewing the assay index
getIntern(tre)  # retrieves data from the @.intern slot
getIntern(tre, "assayIndex")  # get the list item with the associated name from @.intern slot
tre@.intern$assayIndex  # equivalent to above, but without using accessor method

## assays
tre@assays$sensitvity  # the raw format of the assay, indexed by primary key of the same name
tre$sensitivity  # the assay with the row and col data joined to it via the assayIndex

## making a simple summary assay
sens_summary <- tre$sensitivity[,
    .(
        mean_drug1dose=mean(drug1dose),
        mean_drug2dose=mean(drug2dose),
        mean_viability=mean(viability)
    ),
    by=.(drug1id, drug2id, cellid)
]

## assign back to the object
tre$sens_summary <- sens_summary

## view the raw assays
tre@assays$sens_summary  # keyed by sens_summary column
# repeated values in the foreign key of an assay when they are summarized
mutable(tre@.intern$assayIndex)[!is.na(sens_summary), ]  # select rows with sens_summary values
