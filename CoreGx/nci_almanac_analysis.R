library(PharmacoGx)
library(data.table)

# set this to your data directory!
data_dir <- ".local_data"

nci <- readRDS(file.path(data_dir, "NCI_ALMANAC_2017.rds"))

# configure parallelization
nthread <- 4
setDTthreads(nthread)

# -- Fit the Hill curve model to monotherapy drugs and compute associated metrics
bench::system_time({
treatmentResponse(nci) |>
    endoaggregate(
        assay="sensitivity",
        target="mono_raw",
        subset=treatment2id == "",
        viability=mean(viability),
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) |>
    endoaggregate(
        assay="mono_raw",
        target="mono_profiles",
        {
            fit <- PharmacoGx::logLogisticRegression(treatment1dose, viability)
            ic50 <- PharmacoGx::computeIC50(treatment1dose, Hill_fit=fit)
            auc <- PharmacoGx::computeAUC(treatment1dose, Hill_fit=fit,
                area.type="Fitted")
            list(
                HS=fit[['HS']],
                E_inf=fit[['E_inf']] / 100,
                EC50=fit[['EC50']],
                Rsquare=attributes(fit)$Rsquare,
                viability=mean(viability, na.rm=TRUE) / 100,
                auc=auc,
                ic50=ic50
            )
        },
        by=c("treatment1id", "sampleid"),
        enlist=FALSE,
        nthread=nthread
    ) ->
    tre
}) -> m1

# -- Attach the Hill parameters from the monotherapy profiles to the combination data
bench::system_time({
tre |>
    endoaggregate(
        assay="sensitivity",
        target="combo_raw",
        subset=treatment2id != "",
        viability_combo=mean(viability) / 100,
        by=c("treatment1id", "treatment1dose", "treatment2id", "treatment2dose",
            "sampleid")
    ) ->
    tre
}) -> m2

bench::system_time({
tre |>
    mergeAssays(
        "combo_raw",
        "mono_profiles",
        by=c("treatment1id", "sampleid")
    ) |>
    mergeAssays(
        "combo_raw",
        "mono_profiles",
        by.x=c("treatment2id", "sampleid"),
        by.y=c("treatment1id", "sampleid"),
        suffixes=c("_1", "_2")
    ) ->
    tre
}) -> m3

# -- compute our drug synergy metrics
bench::system_time({
tre |>
    endoaggregate(
        assay="combo_raw",
        HSA=min(viability_1, viability_2),
        Bliss=prod(viability_1, viability_2),
        ZIP={
            dose_ratio1 <- (treatment1dose / EC50_1)
            dose_ratio2 <- (treatment2dose / EC50_2)
            ZIP1 <- 1 / (1 + dose_ratio1^HS_1)
            ZIP2 <- 1 / (1 + dose_ratio2^HS_2)
            ZIP1 * ZIP2
        },
        Loewe=PharmacoGx::computeLoewe(
            treatment1dose=treatment1dose,
            treatment2dose=treatment2dose,
            HS_1=HS_1,
            HS_2=HS_2,
            E_inf_1=E_inf_1,
            E_inf_2=E_inf_2,
            EC50_1=EC50_1,
            EC50_2=EC50_2
        ),
        by=assayKeys(tre, "combo_raw")
    ) ->
    ntre1
}) -> m4


# -- compute combination index and score
bench::system_time({
ntre1 |>
    endoaggregate(
        assay="combo_raw",
        HSA_CI=viability_combo / HSA,
        HSA_score=HSA - viability_combo,
        Bliss_CI=viability_combo / Bliss,
        Bliss_score=Bliss - viability_combo,
        ZIP_CI=viability_combo / ZIP,
        ZIP_score=ZIP - viability_combo,
        Loewe_CI=viability_combo / Loewe,
        Loewe_score=Loewe - viability_combo,
        by=assayKeys(ntre1, "combo_raw")
    ) ->
    ntre2
}) -> m5

# -- summarize our combination scores
bench::system_time({
ntre2 |>
    endoaggregate(
        subset=Rsquare_1 > 0.5 & Rsquare_2 > 0.5, # remove bad fits
        assay="combo_raw",
        target="combo_profiles",
        HSA_CI=mean(HSA_CI, na.rm=TRUE),
        Bliss_CI=mean(Bliss_CI, na.rm=TRUE),
        Loewe_CI=mean(Loewe_CI, na.rm=TRUE),
        ZIP_CI=mean(ZIP_CI, na.rm=TRUE),
        HSA_score=mean(HSA_score, na.rm=TRUE),
        Bliss_score=mean(Bliss_score, na.rm=TRUE),
        Loewe_score=mean(Loewe_score, na.rm=TRUE),
        ZIP_score=mean(ZIP_score, na.rm=TRUE),
        by=c("treatment1id", "treatment2id", "sampleid")
    ) ->
    ntre3
}) -> m6

(runtime <- m1 + m2 + m3 + m4 + m5 + m6)

# -- extract our combo scores and have a look
combo_prof <- ntre3$combo_profiles
cor(combo_prof[,
    .(HSA_score, Bliss_score, Loewe_score, ZIP_score)
], method="pearson", use="complete.obs")

# rank by most conservative estimate of synergy
combo_prof[,
    min_score := pmin(HSA_score, Bliss_score, Loewe_score, ZIP_score,
        na.rm=TRUE)
]
combo_prof[,
    max_CI := pmax(HSA_CI, Bliss_CI, Loewe_CI, ZIP_CI, na.rm=TRUE)
]
combo_prof[
    order(-min_score),
    .SD, .SDcols=patterns("*_score|*prop"),
    by=c(assayKeys(ntre3, "combo_profiles"))
] |> head(15)

## _score columns calculate the difference between null model viability and
## observed viability, where positive values are good
combo_prof[,
    lapply(.SD, median, na.rm=TRUE),
    .SDcols=is.numeric,
    by=c("treatment1id", "treatment2id")
][order(-min_score)] |> head(15)


## -- save the results to compare against SYNERGxDB
setorderv(combo_prof, c("min_score"), order=-1L)
fwrite(combo_prof, file=file.path(data_dir, "nci_almanac_synergy_profiles.csv"))