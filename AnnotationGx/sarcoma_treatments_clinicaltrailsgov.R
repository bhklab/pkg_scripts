library(data.table)

sarcoma_therapy_df <- fread("https://clinicaltrials.gov/api/query/study_fields?expr=AREA%5BCondition%5Dsarcoma+AND+AREA%5BIsFDARegulatedDrug%5DYes&fields=InterventionName%2CInterventionOtherName%2CInterventionDescription%2CNCTId%2CBriefTitle%2COfficialTitle%2CCondition&min_rnk=1&max_rnk=1000&fmt=csv")

df_ <- copy(sarcoma_therapy_df)
df_[,
    intervention := lapply(InterventionName,
        FUN=\(x) trimws(strsplit(x, split="\\|")[[1]]))
]
df_[,
    condition := lapply(Condition,
        FUN=\(x) trimws(strsplit(x, split="\\|")[[1]]))
]

flat_df <- df_[,
    lapply(.SD, FUN=unlist),
    .SDcols=c("intervention", "condition", "NCTId", "BriefTitle")
][condition %ilike% "sarcoma", ]

sarc_df <- flat_df[,
    lapply(.SD, \(x) paste0(unique(x), collapse="|")),
    by=NCTId
]