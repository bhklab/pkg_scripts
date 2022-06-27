library(dbplyr)
library(RMariaDB)  # for MySQL DB connector
library(dplyr)
library(data.table)
library(AnnotationGx)

data_dir <- ".local_data"

# =====================
# ---- PharmacoDB query

# -- Connect to the database

# Load credentials from the environment
readRenviron("~/.mysql")
DB_USER <- Sys.getenv("DB_USER")
DB_PASS <- Sys.getenv("DB_PASS")
DB_NAME <- Sys.getenv("DB_NAME")
DB_HOST <- Sys.getenv("DB_HOST")
DB_PORT <- Sys.getenv("DB_PORT")

# Connect to PharmacoDB MySQL Instance
CON <- dbConnect(
    MariaDB(),
    dbname=DB_NAME,
    username=DB_USER,
    password=DB_PASS,
    host=DB_HOST,
    port=DB_PORT
)

# -- Extract the relevant SQL tables
compound <- tbl(CON, "compound")
compound_annot <- tbl(CON, "compound_annotation")

# -- Add the annotations and materialize as a data.table
acompound <- compound |>
    left_join(compound_annot |> select(-id), by=c(id="compound_id")) |>
    as.data.table()
# clean up compound name by removing quotes
acompound[, name := gsub('^\\"*|\\"*$', '', name)]

# close DB conection
dbDisconnect(CON)

# ==========================================================
# ---- Fetch all mechanism of action annotations from ChEMBL
chembl_moa <- getChemblAllMechanisms()

# -- preprocess the table to be writable as a .csv
# extract nested MOA references
mechanism_ref_list <- setNames(
    chembl_moa[["mechanism_refs"]],
    chembl_moa[["record_id"]]
)

# clean up the references to be a single column and aggregate them by record
mechanism_ref_df <- rbindlist(mechanism_ref_list, use.names=TRUE, fill=TRUE,
    idcol="record_id")
mechanism_ref_df[, reference := paste0(ref_type, "=", ref_url)]
mechanism_ref <- mechanism_ref_df[,
    lapply(.SD, function(x) {
        unique_x <- unique(x)
        if (length(unique_x > 1)) {
            list(unique_x)
        } else {
            type.convert(unique_x)
        }
    }),
    by="record_id"
]

chembl_moa[, mechanism_refs := NULL]
chembl_moa_ref <- merge.data.table(
    chembl_moa,
    mechanism_ref[, .(record_id=as.integer(record_id), reference=reference)],
    by="record_id",
    all.x=TRUE
)

# ====================================================================
# ---- Join the mechanisms with PharmacoDB compounds and write to .csv


compound_mechanisms <- merge.data.table(
    acompound,
    chembl_moa_ref,
    by.x="chembl_id",
    by.y="molecule_chembl_id",
    all.x=TRUE
)
compound_mechanisms_2 <- merge.data.table(
    compound_mechanisms[is.na(mechanism_of_action), colnames(acompound),
        with=FALSE],
    chembl_moa_ref,
    by.x="chembl_id",
    by.y="parent_molecule_chembl_id"
)
compounds <- rbind(compound_mechanisms, compound_mechanisms_2, fill=TRUE)
compounds[, reference := unlist(lapply(reference, paste0, collapse="|"))]

# explore common mechanisms of action in PharmacoDB
common_mechanisms <- compounds[
    !is.na(mechanism_of_action),
    .(num_compound=uniqueN(name)),
    by="mechanism_of_action"
][num_compound > 5, ][order(-num_compound), ]

common_mech_compounds <- compounds[
    mechanism_of_action %in% common_mechanisms[["mechanism_of_action"]],
    .(compound=list(unique(name)), ncompound=uniqueN(name)),
    by="mechanism_of_action"
][order(-ncompound)][mechanism_of_action != "Unknown"]

# save common mechanism drugs
fwrite(
    common_mech_compounds,
    file=file.path(data_dir, "common_mechanism_compounds.csv")
)

# save annotated compound table
fwrite(
    compounds,
    file=file.path(data_dir, "pharmacodb_compound_mao.csv")
)
