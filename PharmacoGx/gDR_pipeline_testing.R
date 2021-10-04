#!/bin/R

library(gDR)
library(SummarizedExperiment)
library(BumpyMatrix)
library(qs)
library(data.table)

# 0 -- Load the gDR Style SE
CCLE_gDR <- qread(file.path('local_data', 'CCLE_gDR.qs'))

# 1 -- Average the gDR results
CCLE_averaged <- average_SE(CCLE_gDR, normalized_assay='Viability')

