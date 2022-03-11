library(data.table)
library(qs)
library(MultiAssayExperiment)
library(PharmacoGx)
library(BiocParallel)
library(VennDiagram)
library(UpSetR)
library(ComplexUpset)

## ----- Script parameters

nthread <- 8

# data.table
setDTthreads(nthread)

# BiocParallel
bp <- bpparam()
bpworkers(bp) <- nthread
bpprogressbar(bp) <- TRUE
register(bp)

## ----- Pre-clinical Data

## -- PharmacoGx

sarcsets <- qread("local_data/sarcsets.qs", nthread=nthread)

# microarray

has_microarray <- vapply(sarcsets, \(x) "rna" %in% mDataNames(x), logical(1))
microsets <- sarcsets[has_microarray]

microcells <- lapply(microsets, \(x) colData(molecularProfilesSlot(x)$rna)$cellid)
microsets_sub <- Map(\(x, y) { if (length(y)) subsetTo(x, cells=y, molecular.data.cells=y) else x },
    x=microsets, y=microcells)

microsets_cellInfo <- rbindlist(lapply(microsets_sub, cellInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE
    )[, .SD, .SDcols=patterns("^cellid|^dataset|^cellosarus.*")]

unique_micro_cells_by_disease_dataset_long <- microsets_cellInfo[,
    .(nunique_cells=uniqueN(cellid)),
    by=.(cellosarus.disease, dataset)
][order(-nunique_cells)]

unique_micro_cells_by_disease_dataset <- dcast(
    unique_micro_cells_by_disease_dataset_long,
    cellosarus.disease ~ dataset,
    value.var="nunique_cells",
    fill=0
)
setorderv(unique_micro_cells_by_disease_dataset,
    cols=unique(unique_micro_cells_by_disease_dataset_long$dataset),
    -1L)

fwrite(unique_micro_cells_by_disease_dataset, file=file.path("local_data",
    "sarcsets_unique_micro_cells_by_disease_dataset.csv"))


# rnaseq
has_rnaseq <- vapply(sarcsets, \(x) any(grepl("rnaseq", mDataNames(x))),
    logical(1))
rnaseqsets <- sarcsets[has_rnaseq]

rnaseqcells <- lapply(rnaseqsets,
    \(x) colData(molecularProfilesSlot(x)[[grep("rnaseq", mDataNames(x))[1]]])$cellid)
rnaseqsets_sub <- Map(\(x, y) { if (length(y)) subsetTo(x, cells=y, molecular.data.cells=y) else x },
    x=rnaseqsets, y=rnaseqcells)

rnaseqsets_cellInfo <- rbindlist(lapply(rnaseqsets_sub, cellInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE
    )[, .SD, .SDcols=patterns("^cellid|^dataset|^cellosarus.*")]

unique_rna_cells_by_disease_dataset_long <- rnaseqsets_cellInfo[,
    .(nunique_cells=uniqueN(cellid)),
    by=.(cellosarus.disease, dataset)
][order(-nunique_cells)]

unique_rna_cells_by_disease_dataset <- dcast(
    unique_rna_cells_by_disease_dataset_long,
    cellosarus.disease ~ dataset,
    value.var="nunique_cells",
    fill=0
)
setorderv(unique_rna_cells_by_disease_dataset,
    cols=unique(unique_rna_cells_by_disease_dataset_long$dataset),
    -1L)

fwrite(unique_rna_cells_by_disease_dataset, file=file.path("local_data",
    "sarcsets_unique_rna_cells_by_disease_dataset.csv"))


# number of drugs per dataset by disease
sarcsets_sensInfo <- rbindlist(lapply(sarcsets, sensitivityInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE)[, .(dataset, cellid, drugid)]

sarcsets_cellInfo <- rbindlist(lapply(sarcsets, FUN=cellInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE
    )[, .SD, .SDcols=patterns("^cellid|^dataset|^cellosarus.*")]

sarcsets_cellosaurus_disease <- unique(sarcsets_cellInfo[, .(cellid, cellosarus.disease)])
sarcsets_cells_drugs_disease <- merge.data.table(
    sarcsets_sensInfo, sarcsets_cellosaurus_disease,
    by="cellid",
    all.x=TRUE
    )

sarcsets_drugs_by_disease_long <- sarcsets_cells_drugs_disease[,
    .(ndrugs=uniqueN(drugid)),
    by=.(dataset, cellosarus.disease)
]
sarcsets_drugs_by_disease <- dcast(sarcsets_drugs_by_disease_long,
    cellosarus.disease ~ dataset,
    value.var="ndrugs",
    fill=0
)[!is.na(cellosarus.disease)]
setorderv(sarcsets_drugs_by_disease,
    cols=unique(sarcsets_drugs_by_disease_long$dataset),
    -1L
)

fwrite(sarcsets_drugs_by_disease, file=file.path("local_data",
    "sarcsets_unique_drugs_by_disease_dataset.csv"))






# all molecular data
cellInfo_by_pset <- rbindlist(lapply(sarcsets, FUN=cellInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE
    )[, .SD, .SDcols=patterns("^cellid|dataset|cellosarus.*")]

drugInfo_by_pset <- rbindlist(lapply(sarcsets, FUN=drugInfo), idcol="dataset",
    use.names=TRUE, fill=TRUE)[, .(dataset, drugid, inchikey)]

sarcsets_stats <- cellInfo_by_pset[,
    .(nsamples=length(cellid), nunique_cells=length(unique(cellid)),
        ndatasets=length(unique(dataset))),
    by=cellosarus.disease
][order(-nsamples, -nunique_cells, -ndatasets)]

sarcsets_sample_stats <- cellInfo_by_pset[,
    .(nsamples=length(cellid), nunique_cells=length(unique(cellid)),
        nunique_disease=length(unique(cellosarus.disease))),
    by=dataset
]


datasets_by_disease <- cellInfo_by_pset[,
    .(datasets=paste0(unique(dataset), collapse=";")),
    by=cellosarus.disease
]

disease_by_dataset <- cellInfo_by_pset[,
    .(diseases=paste0(unique(cellosarus.disease), collapse=";")),
    by=dataset
]


# - After Intersecting with TCGA Drugs

sarcsets_tcga_drugs <- qread("local_data/sarcsets_tcga_drugs.qs",
    nthread=nthread)

cellInfo_by_tcga_pset <- rbindlist(lapply(sarcsets_tcga_drugs, FUN=cellInfo),
    idcol="dataset", use.names=TRUE, fill=TRUE
)[, .SD, .SDcols=patterns("^cellid|dataset|cellosarus.*")]

drugInfo_by_pset_tcga <- rbindlist(lapply(sarcsets_tcga_drugs, FUN=drugInfo),
    idcol="dataset", use.names=TRUE, fill=TRUE)[, .(dataset, drugid)]

drugs_by_dataset <- drugInfo_by_pset_tcga[,
    .(ndrugs=length(unique(drugid)), drugs=paste0(unique(drugid), collapse="; ")),
    by=dataset
][order(-ndrugs)]

sarcsets_stats_tcga <- cellInfo_by_tcga_pset[,
    .(nsamples=length(cellid), nunique_cells=length(unique(cellid)),
        ndatasets=length(unique(dataset))),
    by=cellosarus.disease
][order(-nsamples, -nunique_cells, -ndatasets)]

sarcsets_sample_stats_tcga <- cellInfo_by_pset[,
    .(nsamples=length(cellid), nunique_cells=length(unique(cellid)),
        nunique_disease=length(unique(cellosarus.disease))),
    by=dataset
]

datasets_by_disease_tcga <- cellInfo_by_tcga_pset[,
    .(datasets=paste0(unique(dataset), collapse=";")),
    by=cellosarus.disease
]

disease_by_dataset_tcga <- cellInfo_by_tcga_pset[,
    .(diseases=paste0(unique(cellosarus.disease), collapse=";")),
    by=dataset
]

# ----- Clinical Data

gse210 <- qread("local_data/GSE21050_RangedSummarizedExperiment.qs")
gse211 <- qread("local_data/GSE211")
tcga_target <-