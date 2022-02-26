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

# Queries

# -- Cellosaurus and cell tables
cello <- tbl(CON, "cellosaurus")
dataset <- tbl(CON, "dataset")
dataset_cell <- tbl(CON, "dataset_cell")

# -- Build a cellosaurus table with cell and dataset identifiers
cello_df <- cello |>
    left_join(dataset_cell, by=c(cell_id="cell_id")) |>
    left_join(dataset, by=c(dataset_id="id")) |>
    rename(dataset_name="name", cell_name="identifier") |>
    select(-id.x, -id.y) |>
    collect() |>
    as.data.table()

# -- Find all sarcoma cell-lines in all datasets
sarcoma_df <- cello_df |>
    filter(di %like% "sarcoma") |>
    select(cell_name, di, dataset_name)

# -- Group by dataset
sarcoma_cell_df <- sarcoma_df[,
    lapply(.SD, FUN=\(x) list(unique(x))),
    by=cell_name]
sarcoma_cell_df[, di := unlist(di)]
sarcoma_cell_df[, di := gsub("\\|+", "|", di)]

# -- Write table to disk
fwrite(sarcoma_cell_df, file="local_data/pharmacodb_sarcoma_cells.csv")

# -- Filter non-soft-tissue sarcoma
soft_sarcoma_cell_df <- sarcoma_cell_df[!(di %ilike% "osteosarcoma"), ]
fwrite(soft_sarcoma_cell_df,
    file="local_data/pharmacodb_soft_sarcoma_cells.csv")

soft_no_ewing_cell_df <- soft_sarcoma_cell_df[!(di %ilike% "ewing"), ]
fwrite(soft_no_ewing_cell_df,
    file="local_data/pharmacodb_soft_no_ewing_sarcoma_cell.csv")

# -- Mapping cells to drugs
compound <- tbl(CON, "compound")
compound_annot <- tbl(CON, "compound_annotation")
experiment <- tbl(CON, "experiment")
cell <-tbl(CON, "cell")

acompound <- compound |>
    left_join(compound_annot |> select(-id), by=c(id="compound_id")) |>
    select(compound_id=id, name)

compound_cell <- experiment |>
    select(cell_id, compound_id) |>
    distinct() |>
    left_join(cell, by=c(cell_id="id")) |>
    left_join(acompound, by="compound_id") |>
    select(compound=name.y, cell_name=name.x) |>
    collect() |>
    as.data.table()

sarcoma_compound <- merge.data.table(
    sarcoma_df,
    compound_cell,
    by="cell_name"
)

# clean up disease names
sarcoma_compound[, di := gsub(".*; ", "", di)]

sarcoma_compound_by_cell <- sarcoma_compound[,
    .(compounds=list(unique(compound))),
    by=cell_name
]

# load TCGA compound names to intersect
data_dir <- "local_data"
tcga_sarcoma_compound <- fread(file.path(data_dir,
    "tcga_sarcoma_compound_by_sample.csv"))
compound_query <- unique(tcga_sarcoma_compound$drug_name)
pdb_tcga_sarcoma_compound <- sarcoma_compound[compound %in% compound_query, ]

pdb_tcga_compound_by_cell <- pbd_tcga_sarcoma_compound[,
    .(compound=paste(unique(compound), collapse="|")),
    by=cell_name
]

pdb_tcga_compound_by_disease <- pbd_tcga_sarcoma_compound[,
    .(compound=paste(unique(compound), collapse="|"),
    num_cell_line=length(unique(cell_name))),
    by=di
][order(-num_cell_line), ]

fwrite(pdb_tcga_sarcoma_compound, file=file.path(data_dir,
    "pdb_tcga_sarcoma_compound.csv"))
fwrite(pdb_tcga_compound_by_disease, file=file.path(data_dir,
    "pdb_tcga_compound_by_disease.csv"))

sample_compound_disease <- pbd_tcga_sarcoma_compound[di %ilike% "leiomyo",
    .(n_cell_line=length(cell_name),
    n_unique_cell_line=uniqueN(cell_name),
    n_datsets=uniqueN(dataset_name)),
    by=.(di, compound)
][order(-n_cell_line)]

fwrite(sample_compound_disease, file=file.path(data_dir,
    "cell_line_by_compound_disease_pdb.csv"))

sample_compound_lms <- pbd_tcga_sarcoma_compound[di %ilike% "leiomyo",
    .(n_cell_line=length(cell_name),
    n_unique_cell_line=uniqueN(cell_name),
    n_datsets=uniqueN(dataset_name)),
    by=.(compound)
][order(-n_cell_line)]