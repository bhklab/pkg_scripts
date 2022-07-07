library(data.table)

## -- Load and annotate the data

## data can be retrieved from:
## https://codeocean.com/capsule/6322807/tree/v1

## download files:
## 1. data
## - /MAP/NCI-ALMANAC-drugs.txt
## - /MAP/NCI-ALMANAC-samples.txt
## 2. results
## - tbl.nci-almanac.combo_summary.txt
## 3. https://www.synergxdb.ca/dataset_zips/Raw_NCI-ALMANAC.zip
## - Unzip for NCI_compound_names.tsv

# set this to where you downloaded the data
data_dir <- ".local_data"

# actual results
nci <- fread(file.path(data_dir, "tbl.nci-almanac.combo_summary.txt"))
colnames(nci) <- c("sampleid", "treatment1id", "treatment2id", "bliss",
    "loewe", "hsa", "zip")

# annotation files
treatment <- fread(file.path(data_dir, "NCI-ALMANAC_drug.txt"))
sample <- fread(file.path(data_dir, "NCI-ALMANAC_sample.txt"))
nci_treatment <- fread(file.path(data_dir, "NCI_compound_names.tsv"))
names(nci_treatment) <- c("id", "name")
treatment <- merge.data.table(
    treatment,
    nci_treatment[, .(id, treatment_name=name)],
    by.x="nci_id",
    by.y="id"
)[, !"nci_id", with=FALSE]

# add names to treatments

# join to make the results table
nci_synergyxdb <- unique(merge.data.table(
    nci,
    sample[, .(sample_name=standard, id)],
    by.x="sampleid",
    by.y="id",
    allow.cartesian=TRUE
))
nci_synergyxdb <- unique(merge.data.table(
    nci_synergyxdb,
    treatment,
    by.x="treatment1id",
    by.y="synergy_id",
    allow.cartesian=TRUE
))
nci_synergyxdb <- unique(merge.data.table(
    nci_synergyxdb,
    treatment,
    by.x="treatment2id",
    by.y="synergy_id",
    allow.cartesian=TRUE,
    suffixes=c("1", "2")
))
delete_cols <- grep("id", colnames(nci_synergyxdb), value=TRUE)
nci_synergyxdb[, (delete_cols) := NULL]

## -- Save the annotated table to compare with PharmacoGx results
fwrite(nci_synergyxdb, file=file.path(data_dir, "nci_synergxdb.csv"))