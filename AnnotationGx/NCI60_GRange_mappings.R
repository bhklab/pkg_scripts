## ====================
## ---- 6. Using OrgDB

#!/bin/R

library(PharmacoGx)
library(GenomicRanges)
library(GenomicFeatures)
library(SummarizedExperiment)
library(qs)

## =====================
## ---- 0. Load the data
NCI <- qread('NCI60_PSet_Ranges.qs')


# -- 6.1 Load the database package
library(org.Hs.eg.db)
orgDB <- org.Hs.eg.db

for (i in seq_along(molecularProfilesSlot(NCI))) {

    SE <- molecularProfilesSlot(NCI)[[i]]
    rRanges <- rowRanges(SE)

    ## TODO:: Implement mapping for mirna if possible?
    # Skip mirna
    if (metadata(SE)$annotation == 'mirna') next

    # -- 6.2 Try look-up with Symbols
    symbol <- as.character(rRanges$Gene.name)

    # Entrez multimaps to Ensembl gene and trascript, try taking the first result
    cols <- c('GENENAME', 'ENSEMBL')
    res <- as.data.table(select(orgDB, keys=symbol, columns=cols, keytype='SYMBOL'))
    bySymbolDT <- res[, lapply(.SD, first), by=SYMBOL]

    # Join and make sure they match
    rangeMColDT <- as.data.table(mcols(rRanges))
    mergeMColDT <- merge.data.table(rangeMColDT, bySymbolDT,
        by.x='Gene.name', by.y='SYMBOL', all.x=TRUE)

    # Reassign to mcols
    mcols(rRanges) <- as(mergeMColDT, 'DataFrame')

    # -- 6.3 Retry look-up with entrez IDs
    mcolsDT <- as.data.table(mcols(rRanges))
    entrez <- as.character(mcolsDT[is.na(ENSEMBL), Entrez.gene.id])

    res1 <- as.data.table(select(orgDB, keys=entrez, columns=c('ENSEMBL'), 
        keytype='ENTREZID'))
    moreMappingsDT <- res1[!is.na(ENSEMBL), lapply(.SD, first), by=ENTREZID]
    moreMappingsDT[, ENTREZID := as.numeric(ENTREZID)]
    # Do a join with update by reference
    mcolsDT[moreMappingsDT, ENSEMBL := i.ENSEMBL, on=c("Entrez.gene.id==ENTREZID")]

    # Reassign to mcols
    mcols(rRanges) <- as(mcolsDT, 'DataFrame')

    # -- 6.4 Map to ENSEMBL transcripts
    ensembl <- na.omit(mcols(rRanges)$ENSEMBL)

    tx_ids <- as.data.table(select(orgDB, keys=ensembl, columns=c('ENSEMBLTRANS'), 
        keytype='ENSEMBL'))
    txByEnsembl <- tx_ids[!is.na(ENSEMBLTRANS), lapply(.SD, paste, collapse='|'), 
        by=ENSEMBL]

    mcolsDT <- merge.data.table(mcolsDT, txByEnsembl, by='ENSEMBL', all.x=TRUE)
    setnames(mcolsDT,
        old=c('ENSEMBL', 'ENSEMBLTRANS', 'Entrez.gene.id', 'Cytoband', 'Gene.name',
            'Gene.name.url', 'Entrez.gene.id.url', 'Genomic.coordinate.url',
            'GENENAME'),
        new=c('ensembl_gid', 'ensembl_tid', 'entrez_gid', 'cytoband', 'hugo_symbol',
            'gene_name_url', 'entrez_gid_url', 'genomic_coord_url', 
            'gene_description'),
        skip_absent=TRUE)

    mcols(rRanges) <- as(mcolsDT, 'DataFrame')

    # -- 6.5 Assign back to the PSet
    rowRanges(SE) <- rRanges
    molecularProfilesSlot(NCI)[[i]] <- SE

}


##############################
#### PREVIOUS ATTEMPTS #######
##############################

#!/bin/R

library(PharmacoGx)
library(GenomicRanges)
library(GenomicFeatures)
library(SummarizedExperiment)
library(qs)


## =====================
## ---- 0. Load the data
NCI <- qread('NCI60_PSet_Ranges.qs')


## ===============================================
## ---- 1. Extract range information from the PSet
SE <- molecularProfilesSlot(NCI)$rnaseq.iso
rRanges <- rowRanges(SE)


## ==========================
## ---- 2. Using Ensembl TxDB


# -- 2.1 Fetch the Ensembl TxDB object
# 104 is latest version of GRCh37 (hg19)
EnsTxDB <- makeTxDbFromEnsembl(release=104) 
# Match the chromosome naming conventions
seqlevelsStyle(EnsTxDB) <- 'Ensembl'
seqlevelsStyle(rRanges) <- 'Ensembl'

# ---- 2.2 Extract gene and transcript annotations from the TxDB object
tRanges1 <- transcripts(EnsTxDB)
gRanges1 <- genes(EnsTxDB)

# ---- 2.3 Look for overlaps
# To try non-conservative mappings, remove the type argument
olaps_gene1 <- findOverlaps(rRanges, gRanges1, type='within', select='first')
olaps_trans1 <- findOverlaps(rRanges, tRanges1, type='within', select='first')

rRanges1 <- rRanges
mcols(rRanges1) <- cbind(mcols(rRanges1), mcols(gRanges1)[olaps_gene1, ])
mcols(rRanges1) <- cbind(mcols(rRanges1), mcols(tRanges1)[olaps_trans1, ])


## ==========================
## ---- 3. Using Biomart TxDB


# -- 3.1 Fetch the Biomart TxDB object
BmTxDB <- makeTxDbFromBiomart(host="grch37.ensembl.org")
seqlevelsStyle(BmTxDB) <- 'Ensembl'

# -- 3.2 Extract gene and transcript annotations from the TxDB object
tRanges2 <- transcripts(BmTxDB)
gRanges2 <- genes(BmTxDB)

# -- 3.3 Look for overlaps
# olaps_gene2 <- findOverlaps(rRanges, gRanges2, minoverlap=30, select='first')
# olaps_trans2 <- findOverlaps(rRanges, tRanges2, minoverlap=30, select='first')

# Could also try nearest neighbor methods
nolaps_gene2 <- nearest(rRanges, gRanges2)
nolaps_trans2 <- nearest(rRanges, tRanges2)

# -- 3.4 Add overlapping metadata
rRanges2 <- rRanges
mcols(rRanges2) <- cbind(mcols(rRanges2), gene_id=mcols(gRanges2)[nolaps_gene2, ])
mcols(rRanges2) <- cbind(mcols(rRanges2), mcols(tRanges2)[nolaps_trans2, ])


## ==================
## ---- 4. Using UCSC

# -- 4.1 Fetch the TxDB object

# There are a lot of options... which one to choose?
supportedUCSCtables(genome='hg19')

# Trying with UCSC default first
UCSCTxDB <- makeTxDbFromUCSC(genome='hg19', 
    tablename='ensGene')
seqlevelsStyle(UCSCTxDB) <- 'Ensembl'


# -- 4.2 Extract gene and transcript annotations from the TxDB object
tRanges3 <- transcripts(UCSCTxDB)
gRanges3 <- genes(UCSCTxDB)

# -- 4.3 Find the overlaps
olaps_gene3 <- findOverlaps(rRanges, gRanges3, type='within', select='first')
olaps_trans3 <- findOverlaps(rRanges, tRanges3, type='within', select='first')

# -- 4.4 Add the overlap metadata
rRanges3 <- rRanges
mcols(rRanges3) <- cbind(mcols(rRanges3), mcols(gRanges3)[olaps_gene3, ])
mcols(rRanges3) <- cbind(mcols(rRanges3), mcols(tRanges3)[olaps_trans3, ])


## =======================
## ---- 5. Using ensembldb


# v75 corresponds to hg19
library(EnsDb.Hsapiens.v75)
EnsDb <- EnsDb.Hsapiens.v75

# -- 5.1 Extract gene and transcript annotations
tRanges4 <- transcripts(EnsDb)
gRanges4 <- genes(EnsDb)

# -- 5.2 Find overlaps
# olaps_genes4 <- findOverlaps(rRanges, gRanges4, type='within', select='first')
# olaps_trans4 <- findOverlaps(rRanges, tRanges4, type='within', select='first')

# Use nearest neighbor algorithm instead
nolaps_genes4 <- nearest(rRanges, gRanges4)
nolaps_trans4 <- nearest(rRanges, tRanges4)

# -- 5.3 Add the overlapping metadata
rRanges4 <- rRanges
mcols(rRanges4)<- cbind(mcols(rRanges4), mcols(gRanges4)[nolaps_genes4, ])
mcols(rRanges4) <- cbind(mcols(rRanges4), 
    mcols(tRanges4)[nolaps_trans4, c('tx_id', 'tx_biotype')])





