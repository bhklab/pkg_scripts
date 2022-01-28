library(data.table)
library(AnnotationGx)
library(BiocParallel)

synonyms <- fread("local_data/drugid_synonyms.csv")
fda <- getFDAOrangeBookProducts()

# replace special characters with .? to match zero or more single characters
synonyms[, syn := gsub("[^[:alnum:]|]", ".?", syn)]  # exclude | from replacement

bp <- bpparam()
bpprogressbar(bp) <- TRUE
match_by_ingredient <- bplapply(synonyms$unique.drugid,
    FUN=\(drug) unique(fda[
        Ingredient %ilike% synonyms[unique.drugid == drug, ]$syn,
        .(Ingredient, Trade_Name)
    ]),
    BPPARAM=bp
)
names(match_by_ingredient) <- synonyms$unique.drugid
ingredient_df <- rbindlist(match_by_ingredient, idcol="drugid")

match_by_trade_name <- bplapply(synonyms$unique.drugid,
    FUN=\(drug) unique(fda[
        Trade_Name %ilike% synonyms[unique.drugid == drug, ]$syn,
        .(Ingredient, Trade_Name)
    ]),
    BPPARAM=bp
)
names(match_by_trade_name) <- synonyms$unique.drugid
trade_name_df <- rbindlist(match_by_trade_name, idcol="drugid")

# combine the results
fda_df <- setdiff(ingredient_df, trade_name_df)
fda_df <- rbind(ingredient_df, fda_df)
fda_by_drug_df <- fda_df[,
    lapply(.SD, \(x) paste0(unique(x), collapse="|")),
    by=drugid
]
fda_by_drug_df <- merge.data.table(
    fda_by_drug_df,
    synonyms[, .(unique.drugid, syn)],
    by.x="drugid",
    by.y="unique.drugid",
    all.x=TRUE
)

fwrite(fda_by_drug_df, file="local_data/fda_approved_drugids.csv")