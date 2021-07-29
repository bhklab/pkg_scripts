#!/bin/R

library(PharmacoGx)
library(qs)

CCLE <- qread("local_data/CCLE_2015.qs", nthreads=14)

treatmentInfo(CCLE)
treatmentNames(CCLE)

treatments <- treatmentNames(CCLE)[1:5]
CCLE_treat <- subsetByTreatment(CCLE, treatments)

features <- fNames(CCLE, 'rna')[1:5]
CCLE_feat <- subsetByFeature(CCLE, features, 'rna')

samples <- cellNames(CCLE)[1:5]
CCLE_samp <- subsetBySample(CCLE, samples)

LT <- CoreGx:::.sensitivityToLongTable(CCLE)
no_na_assay <- lapply(assays(LT), na.omit)
assays(LT) <- no_na_assay