library(MultiAssayExperiment)
library(httr)
library(jsonlite)
library(data.table)
library(rio)

# -- Fetch ORCESTRA clinical data metadata file
clinical_datasets <- fromJSON(httr::content(GET(
    "https://www.orcestra.ca/api/clinicalgenomics/available"),
    as="text", type="JSON"))
setDT(clinical_datasets)
(clinical_datasets$name)

# -- Extract dataset information and TCGA dataset file
data_dir <- "local_data"
dataset_name <- "TCGA_2.0.1"

tcga_datasets <- clinical_datasets[name == dataset_name, downloadLink[[1]]]
setDT(tcga_datasets)
(tcga_datasets$name)
tcga_name <- "SARC"
file_name <- file.path(data_dir, paste0(dataset_name, "_", tcga_name, ".rds"))
download.file(
    url=tcga_datasets[name == tcga_name, downloadLink],
    destfile=file_name
)

# TODO:: Why does this print weird?
tcga <- readRDS(file_name)

# -- Extract sample metadata
sample_map <- sampleMap(tcga) |>
    as.data.table()
samples_by_assay <- sample_map[,
    .(num_patients=.N, patients=.(primary)),
    by=assay
]
col_data <- colData(tcga) |>
    as.data.table(keep.rownames=TRUE)

all_patients <- col_data[,
    .SD,
    .SDcols=!patterns("\\.")
]

all_assays <- col_data[
    patientID %in% Reduce(f=intersect, samples_by_assay$patients),
    .SD,
    .SDcols=!patterns("\\.")
]

all_patients_by_hist <- all_patients[,
    .(num_patients=.N),
    by=histological_type
]

all_assays_by_hist <- all_assays[,
    .(num_patients=.N),
    by=histological_type
]

for (obj in c("samples_by_assay", "all_patients", "all_assays",
        "all_patients_by_hist", "all_assays_by_hist")) {
    fwrite(
        get(obj),
        file=file.path("local_data", paste0("tcga_sarcoma_", obj, ".csv"))
    )
}

# -- Parsing the colData, which has over 1.4k columns
metadata_df <- data.table(
    colnames(col_data)
)[, tstrsplit(V1, split="\\.")]

## patient.drugs.drug.measure_of_response == RECIST

(unnested_columns <- metadata_df[is.na(V2), unique(V1)])
(nested_once_df <- metadata_df[!is.na(V2) & is.na(V3),
    V2,
    by=V1])
# only patient and admin are nested
(nested_thrice_df <- metadata_df[!is.na(V2) & !is.na(V3) & !is.na(V3) & is.na(V4),
    .(V2, V3), by=V1
])

# admin is has not relevant to treatment response
patient_df <- metadata_df[V1 == "patient", !"V1"]

# within patient, can ignore samples
(patient_nested_once <- patient_df[is.na(V3), unique(V2)])
(patient_nested_twice <- patient_df[!is.na(V3) & is.na(V4), unique(V3), by=V2])

# -- Get TCGA clean drug names from Supplementary Data of Ding et al 2016
# Evaluating the molecule-based prediction of clinical drug responses in cancer
zip_file <- file.path(data_dir, "tcga_response_clean.zip")
download.file(
    "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/19/10.1093_bioinformatics_btw344/4/btw344_supplementary_data.zip?Expires=1648720679&Signature=IKdebruPcYl25vb2JCB6soW5hVUcghXtGsNBa8vsj7SG52EC6RBIRI4Gx~WC60a2NqZOiA~fApeemDzUtGfeEbNXwY6n5FmL-2jg56I8843-L3dPMVKkgh8ArJGohFYQ2rAHi0Hw4-d02i-rRMJzYhI2oPRRZQUf-XrkAYgemLTp4LYjXBypccI5eFYq1EYt2DUa-hFev6LdwwGOK8bu3rmO-NZH-FCTLt7X4h0fFY2z7wQntDAqa~Ki4d3GlW6uRtgczBBQkZANRD6rFWP~K9HNmYFpDBiS~U4j5i1AYqrK1blak-ELsNiSpNX7u9DVLjFYYBv6aEHCU4BrPXEVNg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA",
    destfile=zip_file
)
(unzipped_files <- unzip(zip_file,
    exdir=data_dir
))
tcga_response <- rio::import_list(unzipped_files[2], setclass="data.table")
tcga_drugs <- tcga_response$`Table S1`[-2, ] |>
    setNames(c("tcga_name", "drugbank", "compound_name"))

# -- Remap TCGA names to cleaned up names
tcga_sarcoma_compound <- col_data[!is.na(patient.drugs.drug.drug_name),
    .(sample=rn, drug_name=patient.drugs.drug.drug_name)
]

remapping <- tcga_drugs[!is.na(compound_name) | !is.na(tcga_name),
    paste0(tcga_name, collapse="|"),
    by=compound_name]
for (i in seq_len(nrow(remapping[compound_name != "NA"]))) {
    tcga_sarcoma_compound[drug_name %ilike% remapping[i, ]$V1,
        drug_name := remapping[i, ]$compound_name
    ]
}

# pd 0332991 == Palbociclib, ChEMBL: https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL189963/
tcga_sarcoma_compound[drug_name == "pd 0332991", drug_name := "Palbociclib"]

# Morab-004 == Ontuximab, ChEMBL: https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL2364659/
tcga_sarcoma_compound[drug_name %ilike% "Morab-004", drug_name := "Ontuximab"]

fwrite(tcga_sarcoma_compound,
    file=file.path(data_dir, "tcga_sarcoma_compound_by_sample.csv"))

# Subset TCGA data to only matching drugs for PharmacoDB
compound_tcga_pdb <- fread(file.path(data_dir,
    "pdb_tcga_sarcoma_compound.csv"))

compound_query <- paste(unique(compound_tcga_pdb$compound,
    collapse="|"))
tcga_pdb_sarcoma_compound <- tcga_sarcoma_compound[
    drug_name %in% compound_query,
]

overlap_compound_patients <- tcga_pdb_sarcoma_compound[,
    .(num_patient=length(sample)),
    by=drug_name][order(-num_patient)
]

overlap_compound_disease <- merge.data.table(tcga_pdb_sarcoma_compound,
    all_patients, by.x="sample", by.y="patientID")

patient_compound_disease <- overlap_compound_disease[,
    .(num_patients=length(sample)),
    by=.(drug_name, histological_type)
][order(-num_patients)]

fwrite(patient_compound_disease,
    file=file.path(data_dir, "patient_by_compound_disease_tcga.csv"))