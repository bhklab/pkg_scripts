library(qs)
library(data.table)
library(MultiAssayExperiment)
library(readxl)

# -- Script configuration
nthread <- 8

# --

# load the TCGA.TARGET.GTEx data
mae <- qread("local_data/TCGA.TARGET.GTEx_mae.qs", nthread=nthread)

# load the
target_os_metadata <- readxl::read_excel(
    "metadata/TARGET_OS_ClinicalData_Discovery_20210520.xlsx"
)

# -- Clean up meta data
target_os_meta <- target_os_metadata |>
    setNames(paste0("dimitrios_lab.", gsub(pattern=" ", replacement="_",
        tolower(colnames(target_os_metadata)))))

# -- Join with existing XenaBrowser metadata

# remove tissue sample suffix, so TARGET identifiers match with new metadata
colData(mae[[1]])$patient_id <- gsub("-\\d+$", "", colData(mae[[1]])$sample)
colData(mae[[1]]) <- merge(colData(mae[[1]]), target_os_meta, by.x="patient_id",
    by.y="dimitrios_lab.target_usi", all.x=TRUE)

## ISSUE: there are no Osteosarcoma samples in the XenaTools TARGET data!
## Therefore, we cannot add these treatment annotations!