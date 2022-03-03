library(PharmacoGx)
library(data.table)
library(BiocParallel)

pharmacoset_dir <- file.path("..", "PharmacoGx", "local_data")
if (!dir.exists(pharmacoset_dir)) dir.create(pharmacoset_dir, recursive=TRUE)

pset_metadata <- file.path("..", "AnnotationGx", "local_data")
sarcoma_cells <- fread(file.path(pset_metadata, "pharmacodb_sarcoma_cells.csv")

# clean up disease names
sarcoma_cells[, c("ncit_disease", "ordo_disease") := tstrsplit(di, split="\\|")]
sarcoma_cells[, disease := gsub(".*; ", "", ncit_disease)]

# based on discussion with Dr. Dimitrios Spentos
exclude_diseases <- c("Gliosarcoma", "Pleural sarcomatoid mesothelioma",
    "Uterine carcinosarcoma", "Thyroid gland sarcoma")

has_exclude_diseases <- exclude_diseases %in% sarcoma_cells$disease
if (!all(has_exclude_diseases)) warning("Disease not found: ",
    paste0(exclude_diseases[!has_exclude_diseases], collapse=", "), "!")

sarcoma_df <- sarcoma_cells[!(disease %in% exclude_diseases), ]

# download the relevant PharmacoSets
sarcoma_psets <- sarcoma_df[,
    na.omit(unique(unlist(tstrsplit(dataset_name, split="\\|"))))
]
sarcoma_psets <- gsub("_", ".*", sarcoma_psets)  # fix to regex match psets
avail_psets <- as.data.table(availablePSets(canonical=TRUE))
download_psets <- avail_psets[
    `PSet Name` %ilike% paste0(sarcoma_psets, collapse="|"),
    `PSet Name`
]

bp <- bpparam()
bpworkers(bp) <- 8
bpprogressbar(bp) <- TRUE

sarcoma_psets <- bplapply(download_psets, FUN=downloadPSet,
    saveDir=pharmacoset_dir, BPPARAM=bp)

# subset PharmacoSets to only relevant cell-lines
bp <- bpparam()
bpworkers(bp) <- 8
bpprogressbar(bp) <- TRUE

sarcoma_psets <- bplapply(download_psets, FUN=downloadPSet,
    saveDir=pharmacoset_dir, BPPARAM=bp)

# subset PharmacoSets to only relevant drugs