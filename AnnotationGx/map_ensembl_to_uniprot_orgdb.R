library(org.Hs.eg.db)
library(data.table)

# -- extract our annotations of interest
ensembl <- as.data.table(org.Hs.egENSEMBL)
uniprot <- as.data.table(org.Hs.egUNIPROT)
symbol <- as.data.table(org.Hs.egSYMBOL)
chr_start <- as.data.table(org.Hs.egCHRLOC)
chr_end <- as.data.table(org.Hs.egCHRLOCEND)
ensembl_trans <- as.data.table(org.Hs.egENSEMBLTRANS)

# -- merge genomic coordinations
setnames(chr_end, c("gene_id", "Chromosome"), c("gene_id.end", "Chromosome.end"))
genomic_coord <- cbind(chr_start, chr_end)
# sanity check
stopifnot(genomic_coord[,
    all(gene_id == gene_id.end) & all(Chromosome == Chromosome.end)
])
genomic_coord[, c("gene_id.end", "Chromosome.end") := NULL]

# -- summarize ensembl transcripts
ensembl_trans <- ensembl_trans[,
    .(ensembl_tid=paste0(unique(trans_id), collapse="|")),
    by=gene_id
]

# -- merge with coordinate data
genomic_coord <- merge.data.table(
    genomic_coord,
    ensembl_trans,
    all=TRUE,
    by="gene_id"
)

# -- summarize ensembl transcripts
ensembl <- ensembl[,
    .(ensembl_gid=paste0(unique(ensembl_id), collapse="|")),
    by=gene_id
]

# -- map to genomic coordinates
genomic_coord <- merge.data.table(
    genomic_coord, ensembl, by="gene_id", all=TRUE
)
genomic_coord <- merge.data.table(
    genomic_coord, symbol, by="gene_id", all=TRUE
)

# -- summarize gene uniprot id
uniprot <- uniprot[,
    .(uniprot_id=paste0(unique(uniprot_id), collapse="|")),
    by=gene_id
]

# -- finalize genomic coordinate mapping
genomic_coord <- merge.data.table(
    genomic_coord, uniprot, by="gene_id", all=TRUE
)

# -- save file
fwrite(genomic_coord, file=file.path(".local_data", "org.Hs.eg.db_gene_mappings.csv"))
