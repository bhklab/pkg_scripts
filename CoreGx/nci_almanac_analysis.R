library(PharmacoGx)
library(data.table)
library(BiocParallel)

nci <- readRDS(file.path(".local_data", "NCI_ALMANAC_2017.rds"))

# -- Fit the Hill curve model to monotherapy drugs and compute associated metrics
bench::system_time({
treatmentResponse(nci) |>
    subset(treatment2id == "") |>
    aggregate(
        assay="sensitivity",
        viability=mean(viability),
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) |>
    aggregate2({
        fit <- PharmacoGx::logLogisticRegression(treatment1dose, viability)
        ic50 <- PharmacoGx::computeIC50(treatment1dose, Hill_fit=fit)
        auc <- PharmacoGx::computeAUC(treatment1dose, Hill_fit=fit, area.type="Fitted")
        list(
            HS=fit[['HS']], E_inf=fit[['E_inf']], EC50=fit[['EC50']],
            Rsq=as.numeric(unlist(attributes(fit))),
            auc=auc,
            ic50=ic50
        )},
        by=c("treatment1id", "sampleid"),
        enlist=FALSE,
        nthread=6
    ) -> monotherapy_profiles
})
# Store the results back in our PharmacoSet
treatmentResponse(nci)$monotherapy_profiles <- monotherapy_profiles

# -- Attach the Hill parameters from the monotherapy profiles to the combination data
# Extract the monotherapy fits
treatmentResponse(nci)$monotherapy_profiles |>
    subset(, c("treatment1id", "sampleid", "HS", "E_inf", "EC50")) |>
    unique() ->
    monotherapy_fits

treatmentResponse(nci)$sensitivity |>
    subset(treatment2id != "") ->
    combo_profiles

# Merge with the fits from the monotherapy experiments
combo_profiles <- merge(
    combo_profiles,
    monotherapy_fits,
    by=c("treatment1id", "sampleid")
) # for treatment1id
combo_profiles <- merge(
    combo_profiles,
    monotherapy_fits,
    by.x=c("treatment2id", "sampleid"),
    by.y=c("treatment1id", "sampleid"),
    suffixes=c("_1", "_2")
)
setkeyv(combo_profiles, c("treatment1id", "treatment2id", "sampleid"))
setcolorder(combo_profiles, c("treatment1id", "treatment2id", "sampleid"))

# -- predict viability for each drug in our combination
# fix the scale of E_inf to be a proportion instead of a percent
combo_profiles[,
    c("E_inf_1", "E_inf_2") := .(E_inf_1 / 100, E_inf_2 / 100)
]
# predict viability at the combo doses from the Hill curves fit to monotherapy data
combo_profiles[,
    c("viability_1", "viability_2") := .(
        PharmacoGx:::.Hill(
            log10(treatment1dose),
            c(HS_1, E_inf_1, log10(EC50_1))
        ),
        PharmacoGx:::.Hill(
            log10(treatment1dose),
            c(HS_1, E_inf_2, log10(EC50_2))
        )
    ),
    by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
]

effectToDose <- function(treatment1dose, treatment2dose, viability, 
        E_inf_1, HS_1, EC50_1, E_inf_2, HS_2, EC50_2) {
    (treatment1dose /
        EC50_1 * ((1 - viability) / (viability - E_inf_1))^(1 / HS_1)) +
    (treatment2dose /
        EC50_2 * ((1 - viability) / (viability - E_inf_2))^(1 / HS_2))
}

# -- compute our drug synergy metrics
combo_profiles |>
    aggregate2(
        HSA=min(viability_1, viability_2),
        Bliss=prod(viability_1, viability_2),
        Loewe_CI=(
            (treatment1dose /
                EC50_1 * ((1 - viability) / (viability - E_inf_1))^(1 / HS_1)
            ) +
            (treatment2dose /
                EC50_2 * ((1 - viability) / (viability - E_inf_2))^(1 / HS_2)
            )
        ),
        ZIP_v={
            dose_ratio1 <- (treatment1dose / EC50_1)
            dose_ratio2 <- (treatment2dose / EC50_2)
            ZIP1 <- 1 / (1 + dose_ratio1^HS_1)
            ZIP2 <- 1 / (1 + dose_ratio2^HS_2)
            ZIP1*ZIP2
        },
        ZIP_r={
            dose_ratio1 <- (treatment1dose / EC50_1)
            dose_ratio2 <- (treatment2dose / EC50_2)
            ZIP1 <- dose_ratio1^HS_1 / (1 + dose_ratio1^HS_1)
            ZIP2 <- dose_ratio2^HS_2 / (1 + dose_ratio2^HS_2)
            ZIP1 + ZIP2 - ZIP1*ZIP2
        },
        Loewe=optimise(
                f = obj,
                interval = c(0, 1),
                dose1 = treatment1dose,
                dose2 = treatment2dose,
                HS_1 = HS_1,
                HS_2 = HS_2,
                E_inf_1 = E_inf_1,
                E_inf_2 = E_inf_2,
                EC50_1 = EC50_1,
                EC50_2 = EC50_2
            )$minimum,
        viability_1=viability_1,
        viability_2=viability_2,
        viability=viability / 100,
        nthread=1,
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ) -> combo_profiles1

# -- compute our drug synergy metrics
combo_profiles |>
    aggregate2(
        viability_1=viability_1,
        viability_2=viability_2,
        viability=viability / 100,
        nthread=1,
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ) -> combo_profiles1


# -- compute combination index and score
combo_profiles1 |>
    aggregate2(
        HSA_CI=viability / HSA,
        HSA_score=HSA - viability,
        Bliss_CI=viability / Bliss,
        Bliss_score=Bliss - viability,
        ZIP_CI=viability / ZIP_v,
        ZIP_score=ZIP_v - viability,
        Loewe_score=Loewe - viability,
        Loewe_CI=Loewe_CI,
        viability=viability,
        viability_1=viability_1,
        viability_2=viability_2,
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ) -> combo_profiles2

dt <- merge.data.table(combo_profiles, combo_profiles1, by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid"))


effectToDose <- function(treatment1dose, treatment2dose, viability, 
        E_inf_1, HS_1, EC50_1, E_inf_2, HS_2, EC50_2) {
    (treatment1dose /
        EC50_1 * ((1 - viability) / (viability - E_inf_1))^(1 / HS_1)) +
    (treatment2dose /
        EC50_2 * ((1 - viability) / (viability - E_inf_2))^(1 / HS_2))
}

dt2 <- dt[,
    .(loewe_val=effectToDose(
        treatment1dose, treatment2dose, Loewe,
        E_inf_1, HS_1, EC50_1, 
        E_inf_2, HS_2, EC50_2)), 
    by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")]

combo_profiles2 |>
    aggregate2(lapply(.SD, mean), by=c("treatment1id", "treatment2id", "sampleid"), enlist=FALSE) ->
    combinations

cor(combo_profiles2[, .(HSA_score, Bliss_score, Loewe_score, ZIP_score)], method="spearman")

cor(combinations[, .(HSA_score, Bliss_score, Loewe_score, ZIP_score)], method="spearman")