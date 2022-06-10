library(data.table)
library(checkmate)

## -- Build dummy data for testing

dimensions <- c("treatment", "sample", "dose", "replicate")
dim_data <- lapply(dimensions, FUN=\(x) {
    data_ <- list(LETTERS, rnorm(length(LETTERS))) |>
        setNames(c(paste0(x, "_id"), x))
    dt_ <- as.data.table(data_)
    setkeyv(dt_, paste0(x, "_id"))
    return(dt_)
})
dummy_index <- expand.grid(lapply(dim_data, `[[`, 1)) |>
    setNames(paste0(dimensions, "_id"))
setDT(dummy_index)
dummy_data <- copy(dummy_index)
# Generate some fully connected dummy data
for (i in seq_along(dim_data)) {
    dummy_data <- dummy_data[dim_data[[i]], on=paste0(dimensions[i], "_id")]
}
assay_keys <- compute_assay_keys(dimensions)
split_keys <- lapply(strsplit(assay_keys, "__"), paste0, "_id")
# Add one value for each assay key combination
for (i in seq_along(split_keys)) {
    assay_key <- assay_keys[[i]]
    dummy_data[, c(assay_key) := rnorm(.N), by=c(split_keys[[i]])]
}

# -- Helper function and constructor for star_schema S3 class

compute_assay_keys <- function(dimensions) {
    assay_keys <- unlist(lapply(seq(2, length(dimensions)), FUN=function(n, x) {
        unlist(lapply(combn(x, n, simplify=FALSE), paste0, collapse="__"))
    }, x=dimensions))
    assay_keys
}

find_dimension_metadata <- function(dimensions, data) {
    dim_data <- vector("list", length(dimensions))
    for (i in seq_along(dim_data)) {
        dim_id <- paste0(dimensions[i], "_id")
        column_cardinality <- data[,
            lapply(.SD, uniqueN),
            by=c(dim_id)
        ]
        one_to_one <- which(apply(column_cardinality == 1, 2, FUN=all))
        dim_data[[i]] <- data[, unique(.SD), .SDcols=one_to_one, by=c(dim_id)]
    }
    return(setNames(dim_data, dimensions))
}

extract_fact_data <- function(dimensions, data) {
    assay_keys <- compute_assay_keys(dimensions)
    split_keys <- lapply(strsplit(assay_keys, "__"), paste0, "_id")
    fact_data <- vector("list", length(assay_keys)) |> setNames(assay_keys)
    for (i in seq_along(fact_data)) {
        fact_keys <- split_keys[[i]]
        column_cardinality <- data[,
            lapply(.SD, \(x) uniqueN(x) / .N),
            by=c(fact_keys)
        ]
        one_to_one <- which(apply(column_cardinality == 1, 2, FUN=all))
        fact_data[[i]] <- data[, unique(.SD), .SDcols=one_to_one, by=c(fact_keys)]
    }
    return(setNames(fact_data, assay_keys))
}

star_schema <- function(dimensions, data) {
    dim_data <- find_dimension_metadata(dimensions, data)
    fact_data <- extract_fact_data(dimensions, data)
    structure(list(
        dimensions=dim_data,
        facts=fact_data,
        class="star_schema"
    ))
}