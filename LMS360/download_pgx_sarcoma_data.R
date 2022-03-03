library(PharmacoGx)
library(data.table)
library(BiocParallel)
library(qs)

pharmacoset_dir <- file.path("..", "PharmacoGx", "local_data")
if (!dir.exists(pharmacoset_dir)) dir.create(pharmacoset_dir, recursive=TRUE)

out_dir <- file.path("local_data")

pset_metadata <- file.path("..", "AnnotationGx", "local_data")
sarcoma_cells <- fread(file.path(pset_metadata, "pharmacodb_sarcoma_cells.csv"))

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

pset_file <- file.path(pharmacoset_dir, "sarc_psets.qs")
if (!file.exists(pset_file)) {
    sarc_psets <- bplapply(download_psets, FUN=downloadPSet,
    saveDir=pharmacoset_dir, BPPARAM=bp, timeout=1e6)
    qsave(sarc_psets, file=pset_file, nthread=getDTthreads())
} else {
    sarc_psets <- qread(pset_file, nthread=getDTthreads())
}

# subset PharmacoSets to only relevant cell-lines
keep_cells <- sarcoma_df[, unique(cell_name)]
sarcsets <- bplapply(sarc_psets,
        FUN=\(x, keep_cells) {
    cells <- intersect(cellNames(x), keep_cells)
    subsetTo(x, cell=cells, molecular.data.cells=cells)
}, keep_cells=keep_cells, BPPARAM=bp)

# add additional metadata
colnames(sarcoma_df)[2:ncol(sarcoma_df)] <- paste0("pharmacodb.",
    colnames(sarcoma_df)[2:ncol(sarcoma_df)])
for (i in seq_along(sarcsets)) {
    cellInfo(sarcsets[[i]]) <- merge(cellInfo(sarcsets[[i]]), sarcoma_df,
        by.x="cellid", by.y="cell_name")
}

qsave(sarcsets, file=file.path(out_dir, "sarcsets.qs"), nthread=getDTthreads())

# drugs of interest
