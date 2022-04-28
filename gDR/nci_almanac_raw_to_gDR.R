library(gDRcore)
library(gDR)
library(data.table)

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


nci <- melt.data.table(nci_DT,
    measure.vars=c("TESTVALUE", "CONTROLVALUE"),
    variable.name="condition",
    value.name="optical_density"
)

# add duration data
nci[, duration := 2]

# format treated vs untreated label
nci[condition == "TESTVALUE", condition := "treated"]
nci[condition != "treated", condition := "untreated"]

# configure column names to match gDR identifiers
gdr_names <- c(
    cellid=get_env_identifiers("cellline"),
    drug1dose=get_env_identifiers("concentration"),
    drug2dose=get_env_identifiers("concentration2"),
    drug1id=get_env_identifiers("drug_name"),
    drug2id=get_env_identifiers("drug_name2"),
    NSC1=get_env_identifiers("drug"),
    NSC2=get_env_identifiers("drug2"),
    PANEL=get_env_identifiers("cellline_tissue"),
    optical_density="ReadoutValue",
    PLATE=get_env_identifiers("barcode")[2],
    duration=get_env_identifiers("duration"),
    TZVALUE="BackgroundValue"
)

setnames(nci, names(gdr_names), gdr_names)

nci_small <- nci[
    clid %in% unique(clid)[1:5] & 
        DrugName %in% unique(DrugName)[1:5] & 
        DrugName_2 %in% unique(DrugName_2)[1:5],
]

# unmatch controls to see if the create_SE function works
nci_small[condition == "untreated", c("DrugName", "DrugName_2") := "untreated"]

se <- create_SE(nci_small, nested_identifiers="Plate")
