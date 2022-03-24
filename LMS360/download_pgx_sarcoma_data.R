library(PharmacoGx)
library(data.table)
library(BiocParallel)
library(qs)
library(AnnotationGx)  # remotes::install_github("bhklab/AnnotationGx")

# -- configure script paths and parameters

# paths
pharmacoset_dir <- file.path("..", "PharmacoGx", "local_data")
if (!dir.exists(pharmacoset_dir)) dir.create(pharmacoset_dir, recursive=TRUE)
out_dir <- file.path("local_data")

# parameters
nthread <- 8

# -- load PharmacoDB cell-line metadata
pset_metadata <- file.path("..", "AnnotationGx", "local_data")
sarcoma_cells <- fread(file.path(pset_metadata, "pharmacodb_sarcoma_cells.csv"))

# clean up disease names
sarcoma_cells[, c("ncit_disease", "ordo_disease") := tstrsplit(di, split="\\|")]
sarcoma_cells[, disease := gsub(".*; ", "", ncit_disease)]

# remove mislabelled sarcomas based on discussion with Dr. Dimitrios Spentos
exclude_diseases <- c("Gliosarcoma", "Pleural sarcomatoid mesothelioma",
    "Uterine carcinosarcoma", "Thyroid gland sarcoma",
    "Ovarian clear cell adenocarcinoma", "Rhabdoid tumor of the kidney",
    "Lung adenocarcinoma", NA_character_)

has_exclude_diseases <- exclude_diseases %in% sarcoma_cells$disease
if (!all(has_exclude_diseases)) warning("Disease not found: ",
    paste0(exclude_diseases[!has_exclude_diseases], collapse=", "), "!")

sarcoma_df <- sarcoma_cells[!(disease %in% exclude_diseases), ]

# -- download the relevant PharmacoSets
sarcoma_psets <- sarcoma_df[!is.na(disease),
    na.omit(unique(unlist(tstrsplit(dataset_name, split="\\|"))))
]
sarcoma_psets <- gsub("_", ".*", sarcoma_psets)  # fix to regex match psets
avail_psets <- as.data.table(availablePSets(canonical=TRUE))
download_psets <- avail_psets[
    `PSet Name` %ilike% paste0(sarcoma_psets, collapse="|"),
    `PSet Name`
]

# configure parallelization settings for BiocParallel::bplapply
bp <- bpparam()
bpworkers(bp) <- nthread
bpprogressbar(bp) <- TRUE

pset_file <- file.path(pharmacoset_dir, "sarc_psets.qs")
if (!file.exists(pset_file)) {
    sarc_psets <- bplapply(download_psets, FUN=downloadPSet,
    saveDir=pharmacoset_dir, BPPARAM=bp, timeout=1e6)
    qsave(sarc_psets, file=pset_file, nthread=nthread)
} else {
    sarc_psets <- qread(pset_file, nthread=nthread)
}

# -- subset PharmacoSets to only relevant cell-lines
keep_cells <- sarcoma_df[, unique(cell_name)]
sarcsets <- bplapply(sarc_psets,
        FUN=\(x, keep_cells) {
    cells <- intersect(cellNames(x), keep_cells)
    subsetTo(x, cell=cells, molecular.data.cells=cells)
}, keep_cells=keep_cells, BPPARAM=bp)
sarcsets <- setNames(sarcsets, download_psets)

# -- Add NCI Sarcoma PSet to list
# FIXME:: remove this when NCI Sarcoma gets in ORCESTRA
nci <- readRDS(file.path(out_dir, "NCI_Sarcoma.rds"))

sarcsets <- c(sarcsets, list(NCI_Sarcoma=nci))

# -- add additional Cellosarus metadata
cells <- unique(Reduce(c, lapply(sarcsets, FUN=cellNames)))
cello_df <- getCellosaurus(cells)  # downloads from Cellosaurus

# add additional PharmacoDB metadata
colnames(cello_df)[2:ncol(cello_df)] <- paste0("cellosaurus.",
    colnames(cello_df)[2:ncol(cello_df)])  # label columns from PharmacoDB
for (i in seq_along(sarcsets)) {
    cellInfo(sarcsets[[i]]) <- merge(cellInfo(sarcsets[[i]]), cello_df,
        by.x="cellid", by.y="standard_name", all.x=TRUE)
    # fix rownames dropped by join
    rownames(cellInfo(sarcsets[[i]])) <- cellInfo(sarcsets[[i]])$cellid
}

# fix molecular data
molecularProfilesSlot(sarcsets$NCI_Sarcoma)$rna <- as(
    molecularProfilesSlot(sarcsets$NCI_Sarcoma)$rna,
    "SummarizedExperiment"
)

# remove diseases flagged as not sarcoma
cInfo <- as.data.table(cellInfo(sarcsets[["NCI_Sarcoma"]]))
keep_cells <- cInfo[
    !(cellosaurus.disease %in% exclude_diseases | is.na(cellosaurus.disease)),
    unique(cellid)
]
sarcsets$NCI_Sarcoma <- subsetTo(sarcsets$NCI_Sarcoma, cells=keep_cells,
    molecular.data.cells=keep_cells)

qsave(sarcsets, file=file.path(out_dir, "sarcsets.qs"),
    nthread=nthread)

# -- subset PharmacoSets to only drugs also in TCGA
tcga_pgx_df <- fread(file.path(pset_metadata,
    "pdb_tcga_compound_by_disease.csv"))
tcga_pgx_drugs <- tcga_pgx_df[,
    na.omit(unique(unlist(tstrsplit(compound, split="\\|"))))
]
# remove drugs with weird regex match
tcga_pgx_drugs <- tcga_pgx_drugs[1:9]
fwrite(list(compound=tcga_pgx_drugs), file=file.path(out_dir,
    "tcga_pgx_drugs.csv"))

sarcsets_tcga_drugs <- bplapply(sarcsets,
        FUN=\(x, keep_drugs) {
    drugs <- intersect(drugNames(x), keep_drugs)
    subsetTo(x, drugs=drugs)
}, keep_drugs=tcga_pgx_drugs, BPPARAM=bp)

qsave(sarcsets_tcga_drugs, file=file.path(out_dir, "sarcsets_tcga_drugs.qs"),
    nthread=nthread)