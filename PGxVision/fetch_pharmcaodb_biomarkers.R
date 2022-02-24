library(dbplyr)
library(RMariaDB)  # for MySQL DB connector
library(dplyr)
library(data.table)

# Load credentials from the environment
readRenviron("~/.mysql")
MYSQL_USER <- Sys.getenv("MYSQL_USER")
MYSQL_PW <- Sys.getenv("MYSQL_PW")
HOST <- Sys.getenv("PHARMACODB")
PORT <- Sys.getenv("PHARMACODB_PORT")

# Connect to PharmacoDB MySQL Instance
CON <- dbConnect(
    MariaDB(),
    dbname="pharmacodb_test",
    username=MYSQL_USER,
    password=MYSQL_PW,
    host=HOST,
    port=PORT
)

# Get tables of interest
gene <- tbl(CON, "gene")
gene_annot <- tbl(CON, "gene_annotation")
compound <- tbl(CON, "compound")
compound_annot <- tbl(CON, "compound_annotation")
tissue <- tbl(CON, "tissue")
biomarker <- tbl(CON, "gene_compound_tissue")

# Join
tissue_df <- tissue |>
    rename(tissue="name")

gene_df <- gene_annot |>
    filter(symbol != "NA") |>
    select(gene_id, symbol) |>
    left_join(gene, by=c(gene_id="id")) |>
    rename(gene_symbol="symbol", ensembl_id="name")


compound_df <- compound_annot |>
    select(compound_id, inchikey, smiles, pubchem, chembl_id) |>
    left_join(compound, by=c(compound_id="id")) |>
    select(-compound_uid) |>
    rename(compound_name="name")

biomarker_df <- biomarker |>
    filter(pvalue < 0.05) |>
    select(gene_id, compound_id, tissue_id, estimate, pvalue) |>
    inner_join(gene_df, by="gene_id") |>
    inner_join(compound_df, by="compound_id") |>
    inner_join(tissue, by=c("tissue_id"="id")) |>
    select(-gene_id, -compound_id, -tissue_id) |>
    collect()

fwrite(biomarker_df, file=file.path("local_data",
    "pharmacodb_biomarkers.csv"))