library(PharmacoGx)
library(data.table)

pset_name <- "CCLE_2015"
data_dir <- ".local_data"
pset_path <- file.path(data_dir, paste0(pset_name, ".rds"))

if (!file.exists(pset_path)) {
    pset <- downloadPSet(pset_name, saveDir=dirname(pset_path))
} else {
    pset <- readRDS(pset_path)
}

sr <- copy(sensitivityRaw(pset))
si <- copy(sensitivityInfo(pset))
sp <- copy(sensitivityProfiles(pset))

# fix issue where some rownames mismatch
mismatch_rownames <- rownames(si) != rownames(sr)
mismatch_rownames <- mismatch_rownames & (rownames(sp) != rownames(si))
if (any(mismatch_rownames)) {
    warning("Mismatching rownames! Fixing pset...")
    si <- si[!mismatch_rownames, ]
    sr <- sr[!mismatch_rownames, , ]
    sp <- sp[!mismatch_rownames, ]
    sensitivityInfo(pset) <- si
    sensitivityRaw(pset) <- sr
    sensitivityProfiles(pset) <- sp
}

CoreGx:::.compareSensitivityInfo(pset)

srDT <- as.data.table(sr, keep.rownames="rownames")