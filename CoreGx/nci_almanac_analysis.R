library(PharmacoGx)
library(data.table)
library(BiocParallel)

nci <- readRDS(file.path(".local_data", "NCI_ALMANAC_2017.rds"))


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
                auc=auc,
                ic50=ic50
            )
        },
        by=c("treatment1id", "sampleid"),
        enlist=FALSE,
        nthread=22
    ) ->
    tre
}) -> m1

bench::system_time({
# -- Attach the Hill parameters from the monotherapy profiles to the combination data
# Extract the monotherapy fits
tre |>
    endoaggregate(
        assay="sensitivity",
        target="combo_raw",
        subset=treatment2id != "",
        viability=mean(viability) / 100,
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
        by=c(treatment2id="treatment1id", "sampleid"),
        suffixes=c("_1", "_2")
    ) ->
    tre
}) -> m3

# -- predict viability for each drug in our combination

bench::system_time({
tre |>
    endoaggregate(
        assay="combo_raw",
        target="combo_raw",
        {
            v1 <- PharmacoGx:::.Hill(log10(treatment1dose),
                c(HS_1, E_inf_1, log10(EC50_1)))
            v2 <-  PharmacoGx:::.Hill(log10(treatment2dose),
                c(HS_2, E_inf_2, log10(EC50_2)))
            list(viability_1=v1, viability_2=v2)
        },
        by=assayKeys(tre, "combo_raw"),
        enlist=FALSE
    ) ->
    ntre1
}) -> m4

# -- compute our drug synergy metrics
bench::system_time({
ntre1 |>
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
        by=assayKeys(ntre1, "combo_raw")
    ) ->
    ntre2
}) -> m5


# -- compute combination index and score
bench::system_time({
ntre2 |>
    endoaggregate(
        assay="combo_raw",
        HSA_CI=viability / HSA,
        HSA_score=HSA - viability,
        Bliss_CI=viability / Bliss,
        Bliss_score=Bliss - viability,
        ZIP_CI=viability / ZIP,
        ZIP_score=ZIP - viability,
        Loewe_score=Loewe - viability,
        Loewe_CI=viability / Loewe,
        by=assayKeys(ntre2, "combo_raw")
    ) ->
    ntre2
}) -> m6

# -- summarize our combination scores
bench::system_time({
ntre2 |>
    endoaggregate(
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
}) -> m7

# -- extract our combo scores and have a look
combo_prof <- ntre3$combo_profiles
cor(combo_prof[,
    .(HSA_score, Bliss_score, Loewe_score, ZIP_score)
], method="pearson", use="complete.obs")

## _score columns calculate the difference between null model viability and
## observed viability, where positive values are good
combo_prof[
    order(-Bliss_score),
    lapply(.SD, mean, na.rm=TRUE),
    .SDcols=is.numeric,
    by=c("treatment1id", "treatment2id")
]
