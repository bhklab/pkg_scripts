library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(data.table)

# -- Fetch metadata
data(diseaseCodes, package="TCGAutils")


# -- Download all Sarcoma data into an MAE
sarc_mae <- curatedTCGAData(diseaseCode="SARC",
    c("RNASeq2Gene", "Methylation", "*CNA*", "*miRNASeq*", "*RPPA*"),
    version="2.0.1",
    dry.run=FALSE
)

sample_map <- sampleMap(sarc_mae) |>
    as.data.table()
samples_by_assay <- sample_map[,
    .(num_patients=.N, patients=.(primary)),
    by=assay
]
col_data <- colData(sarc_mae) |>
    as.data.table()

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