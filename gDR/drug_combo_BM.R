#!/bin/R

# -- Required
library(BumpyMatrix)
library(gDR)
library(SummarizedExperiment)

# ---- 
testFiles <- list.files('inst/testdata', full.names=TRUE)

SE <- readRDS(testFiles[6])