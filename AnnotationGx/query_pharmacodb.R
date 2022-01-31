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

# -- Cellosarus table
cello <- tbl(CON, "cellosaurus")
cell <- tbl(CON, "cell")

cell_df <- cello |>
    left_join(cell, by=c(cell_id="id")) |>
    filter(di %like% "%sarcoma%") |>
    select(name, di) |>
    collect()