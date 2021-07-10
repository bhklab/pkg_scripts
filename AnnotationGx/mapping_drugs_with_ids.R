#!/bin/R

library(AnnotationGx)
library(data.table)


# ---- Annotating lab standardized compound .csv

# -- Load data
compound <- fread('local_data/drugs_with_ids.csv')
ids <- unique(na.omit(compound$unique.drugid))
inchikeys <- unique(na.omit(compound$inchikey))
smiles <- unique(na.omit(compound$smiles))

# -- SMILES -> CID
#> This is really slow
t1b <- Sys.time()
CIDfromSMILES <- getPubChemCompound(smiles[1:5], from='fastidentity/smiles', 
    to='cids', proxy=TRUE, batch=FALSE)
t2b <- Sys.time()
q2_time <- t2b - t1b
q2_time

# -- CID -> properties
cids <- na.omit(unique(CIDfromSMILES$cids))
NameFromCID <- getPubChemCompound(cids[1:100], proxy=TRUE)


# ---- Other useful cases

# -- inchikeys -> CID
# Note: many CID per inchikey, not that useful
t1a <- Sys.time()
CIDfromInchikey <- getPubChemCompound(inchikeys[1:100], from='inchikey', to='cids', 
    proxy=TRUE, batch=FALSE)
t2a <- Sys.time()
q1_time = t2a - t1a
q1_time

# -- matching syonyms
SynonymsAndIDs <- getPubChemAnnotations('Synonyms and Identifiers', proxy=TRUE)
`%vlike%` <- function(patterns, vector) which(vapply(patterns, FUN=like, vector=vector, logical(length(vector))))
matchingSynoyms <- SynonymsAndIDs[, .(matches=which(compound$unique.drugid %ilike% Synonyms)), by=CID]
