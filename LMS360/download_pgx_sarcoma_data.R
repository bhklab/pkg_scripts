library(PharmacoGx)
library(data.table)

pharmacoset_dir <- file.path("..", "PharmacoGx", "local_data")
if (!dir.exists(pharmacoset_dir)) dir.create(pharmacoset_dir, recursive=TRUE)

pset_metadata <- file.path("..", "AnnotationGx", "local_data")
sarcoma_cells <- fread(file.path(pset_metadata, "pharmacodb_sarcoma_cells.csv")

# clean up disease names
sarcoma_cells[, disease := gsub(".*; ", "", di)]

# based on discussion with Dr. Dimitrios Spentos
exclude_diseases <-