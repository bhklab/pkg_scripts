library(data.table)

drugs <- fread("https://github.com/BHKLAB-Pachyderm/Annotations/raw/master/drugs_with_ids.csv")

.collapse <- function(...) paste0(..., sep="|", collapse="|")
syn <- melt(drugs,
    id.vars="unique.drugid",
    measure.vars=patterns(".*drugid$"),
    na.rm=TRUE
)[, !"variable"]
synonyms <- syn[!(value %like% "///" | value == ""),
    .(syn=paste0(unique(na.omit(value)), collapse="|")),
    by=unique.drugid
]
synonyms[, syn := gsub("\\|\\|", "|", syn)]

fwrite(synonyms, file="local_data/drugid_synonyms.csv")