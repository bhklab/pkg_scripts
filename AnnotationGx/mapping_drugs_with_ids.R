#!/bin/R

library(AnnotationGx)
library(data.table)
library(BiocParallel)
library(doParallel)
library(bench)
library(async)
library(qs)

## ==============================================
## ---- Annotating lab standardized compound .csv

# -- Load data
compound <- fread('local_data/drugs_with_ids.csv')
ids <- unique(na.omit(compound$unique.drugid))
inchikeys <- unique(na.omit(compound$inchikey))
smiles <- unique(na.omit(compound$smiles))

## 1. Map from names via PubChem Compound API
CIDfromNames <- getPubChemCompound(ids[1:1000], from='name', 
    to='cids', proxy=TRUE)


# find the last batch number written to disk
offset <-  max(as.numeric(gsub('^CIDfromSMILES_|.qs$', '', 
        list.files('local_data', pattern='CIDfromSMILES_')))) + 1

# -- SMILES -> CID
# Note: This is really slow
t1b <- Sys.time()
CIDfromSMILES <- getPubChemCompound(smiles[1:100], from='fastidentity/smiles', 
    to='cids', proxy=TRUE, batch=FALSE)
t2b <- Sys.time()
q2_time <- t2b - t1b
q2_time

# -- CID -> properties
cids <- na.omit(unique(CIDfromSMILES$cids))
NameFromCID <- getPubChemCompound(cids, proxy=TRUE, 
    properties=c("Title", "IUPACName"))

# -- Load the data an join together with drugs_with_ids
NameFromCID <- fread('local_data/NameFromSMILES.csv')
setnames(NameFromCID, c('cids', 'Title'), c('new_cid', 'PubChem_title'))
drugs_with_ids <- fread('local_data/drugs_with_ids.csv')
drugs_with_new_ids <- merge.data.table(drugs_with_ids, NameFromCID, by='smiles', 
    all.x=TRUE)

# Check for different CIDs
drugs_with_new_ids[!is.na(new_cid), sum(cid != new_cid, na.rm=TRUE)]

# -- Get an identifier for all the names we couldn't map
inchikeys <- drugs_with_new_ids[is.na(new_cid) & !is.na(inchikey), inchikey]
drug_names <- drugs_with_new_ids[is.na(new_cid) & is.na(inchikey), unique.drugid]

# -- Map missing new_cid from inchikey
CIDfromInchikey <- getPubChemCompound(inchikeys, from='inchikey', batch=FALSE, 
    proxy=TRUE)


# t1b <- Sys.time()
# CIDfromSMILES <- getPubChemCompound(smiles[1:100], from='fastidentity/smiles', 
#     to='cids', proxy=TRUE, batch=FALSE, BPPARAM=bp)
# t2b <- Sys.time()
# q2_time <- t2b - t1b
# q2_time

# -- Testing out async package for faster queries

# proxyDT <- AnnotationGx:::proxyManager$get_proxies()
# proxyDT <- proxyDT[, url := paste0(ip, ":", port), by=.(ip, port)]
# proxy_urls <- proxyDT$url

# queries <- unlist(getPubChemCompound(smiles[1:1000], from='fastidentity/smiles', 
#     to='cids', proxy=TRUE, batch=FALSE, query_only=TRUE, BPPARAM=bp))
# batch_size <- ceiling(length(queries) / nthreads)
# split_queries <- split(queries, floor(seq_along(queries) / batch_size))

#' @description
#' An async function that will retry ten random proxies until a result is
#' retrieved
#'
#' @return `deferred` R6 object. Call `synchronize` to evaluate the result
#'   in an async context.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom async async_retry http_get
# http_get_proxy_retry <- async(function(url, proxies) {
#     async_retry(
#         function() http_get(url, options=list(proxy=sample(proxies, 1))),
#         times=5
#     )$then(function(x) jsonlite::fromJSON(rawToChar(x$content)))$
#       catch(error=function(e) { print(e); cat('\n'); return(url)})
# })

# batch_time <- bench_time({
#     batch_result <- synchronise(async_map(queries[1:100], 
#         .f=http_get_proxy_retry, proxies=proxy_urls))
# })




# ---- Other useful cases

# -- inchikeys -> CID
# Note: many CID per inchikey, not that useful
t1a <- Sys.time()
CIDfromInchikey <- getPubChemCompound(inchikeys[1:1000], from='inchikey', to='cids', 
    proxy=TRUE, batch=FALSE)
t2a <- Sys.time()
q1_time = t2a - t1a
q1_time

# -- matching syonyms
# Note: Not working correctly
SynonymsAndIDs <- getPubChemAnnotations('Synonyms and Identifiers', proxy=TRUE)
`%vlike%` <- function(patterns, vector) which(vapply(patterns, FUN=like, vector=vector, logical(length(vector))))
# matchingSynoyms <- SynonymsAndIDs[, .(matches=which(compound$unique.drugid %ilike% Synonyms)), by=CID]
