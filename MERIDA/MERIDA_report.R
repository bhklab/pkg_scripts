library(data.table)
library(ggplot2)

# -- compile the model results
eval_files <- list.files("results", pattern="MERIDA.*model_eval_.*.csv",
    full.names=TRUE)
eval_data <- lapply(eval_files, FUN=fread)
model_eval <- rbindlist(setNames(eval_data, eval_files),
    idcol="src_file")
# clean up file paths to extract dataset and drug names
model_eval[,
    src_file := gsub("^results/MERIDA_|_model_eval_t.*csv$", "", src_file)
]
model_eval[,
    c("dataset", "drug") := tstrsplit(src_file, "_")
]
model_eval[, training_set := FALSE]
model_eval[dataset == "CCLE", training_set := TRUE]

# -- subset to relevant metrics and clean up for plotting
df_ <- model_eval[
    fold == "all",
    .(dataset, drug, sensitivity, resistance, mcc, ppv=`Pos Pred Value`,
        npv=`Neg Pred Value`, bal_acc=`Balanced Accuracy`)
]
plot_df <- df_[, .SD, .SDcols=!c("sensitivity", "resistance")]
plot_df_melt <- melt(plot_df,
    id.vars=c("dataset", "drug"),
    variable.name="metric", value.name="value"
)
plot_df_melt[,
    metric_pretty := factor(fcase(
        metric == "mcc", "Matthews Correlation Coefficient",
        metric == "ppv", "Positivie Predictive Value",
        metric == "npv", "Negative Predictive Value",
        metric == "bal_acc", "Balanced Accuracy"
    ))
]
plot_df_melt[, step := fifelse(dataset == "CCLE", "Training", "Testing")]

df_ <- plot_df_melt[
    metric == "mcc",
    if (min(value) > 0.2) .SD, ],
    by=drug
]
# -- plot
grouped_bar_plot <- ggplot(df_, aes(x=drug, y=value, fill=dataset)) +
    geom_bar(position="dodge", stat="identity") +
    coord_flip()
