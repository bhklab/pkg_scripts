library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)
library(qs)

## NOTE:: This script requires ~ 40 GB of RAM

# -- Find datasets of interest

# Load dataset metadata index
data(XenaData)
xd <- as.data.table(XenaData)

(xena_hosts <- unique(xd$XenaHostNames))
host <- "publicHub"

(target_cohorts <- unique(xd[XenaCohorts %ilike% "^TARGET", ]))
