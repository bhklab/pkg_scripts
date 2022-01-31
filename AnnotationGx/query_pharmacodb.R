library(dbplyr)
library(RMariaDB)  # for MySQL DB connector
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
