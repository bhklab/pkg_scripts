#!/bin/R

library(AnnotationGx)
library(data.table)

# ---- Annotating lab standardized compound .csv
compound <- fread('local_data/drugs_with_ids.csv')
ids <- compound$unique.drugid
inchikeys <- compound$inchikey
CIDfromInchikey <- getPubChemCompound(inchikeys, from='inchikey', to='cids', 
    proxy=TRUE)
