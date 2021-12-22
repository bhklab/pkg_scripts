#!/bin/R

library(gDR)
library(SummarizedExperiment)
library(BumpyMatrix)
library(qs)
library(data.table)

# 0 -- Load the gDR Style SE
CCLE_SE <- qread(file.path('local_data', 'CCLE_gDR.qs'))


# 1 -- Customize the identifiers to match our data
pgx_identifiers <- c(
    drug="drug1id",
    concentration="drug1dose",
    cellline="cellid"
)
for (i in seq_along(pgx_identifiers)) {
    gDRutils::set_env_identifier(
        names(pgx_identifiers)[i],
        pgx_identifiers[i]
    )
}

# 1 -- Average the gDR results
CCLE_averaged <- average_SE(CCLE_SE)
