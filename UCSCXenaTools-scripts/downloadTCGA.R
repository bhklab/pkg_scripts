library(UCSCXenaTools)
library(dplyr)
library(data.table)
library(MultiAssayExperiment)
library(S4Vectors)

# -- Find datasets of interest

# Load dataset metadata index
data(XenaData)
xd <- as.data.table(XenaData)

(xena_hosts <- unique(xd$XenaHostNames))


(tcga_cohorts <- unique(xd[XenaCohorts %ilike% "^TCGA", ]$XenaCohorts))