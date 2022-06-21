library(PharmacoGx)
library(data.table)
library(BiocParallel)

nci <- readRDS(file.path(".local_data", "NCI_ALMANAC_2017.rds"))
.nci <- copy(nci)

# -- Fit the Hill curve model to monotherapy drugs and compute associated metrics
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
            Rsquare=attributes(fit)$Rsquare,
            auc=auc,
            ic50=ic50
        )},
        by=c("treatment1id", "sampleid"),
        enlist=FALSE,
        nthread=20
    ) ->
    monotherapy_profiles
# Store the results back in our PharmacoSet
treatmentResponse(nci)$monotherapy_profiles <- monotherapy_profiles

# -- Attach the Hill parameters from the monotherapy profiles to the combination data
# Extract the monotherapy fits
treatmentResponse(nci)$monotherapy_profiles |>
    subset(Rsquare > 0.5, c("treatment1id", "sampleid", "HS", "E_inf", "EC50")) |>
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

# -- compute our drug synergy metrics
combo_profiles |>
    aggregate2(
        HSA=min(viability_1, viability_2),
        Bliss=prod(viability_1, viability_2),
        Loewe_CI=(
            (treatment1dose /
                EC50_1 * ((1 - viability_1) / (viability_1 - E_inf_1))^(1 / HS_1)
            ) +
            (treatment2dose /
                EC50_2 * ((1 - viability_2) / (viability_2 - E_inf_2))^(1 / HS_2)
            )
        ),
        ZIP={
            dose_ratio1 <- (treatment1dose / EC50_1)
            dose_ratio2 <- (treatment2dose / EC50_2)
            ZIP1 <- 1 / (1 + dose_ratio1^HS_1)
            ZIP2 <- 1 / (1 + dose_ratio2^HS_2)
            ZIP1 * ZIP2
        },
        ZIP2={
            dose_ratio1 <- (treatment1dose / EC50_1)
            dose_ratio2 <- (treatment2dose / EC50_2)
            ZIP1 <- E_inf_1 * dose_ratio1^HS_1 / (1 + dose_ratio1^HS_1)
            ZIP2 <- E_inf_2 * dose_ratio2^HS_2 / (1 + dose_ratio2^HS_2)
            ZIP1 * ZIP2
        },
        PANEL=PANEL,
        viability_1=viability_1,
        viability_2=viability_2,
        viability=viability / 100,
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ) ->
    combo_profiles1

# -- compute combination index and score
combo_profiles1 |>
    aggregate2(
        HSA_CI=viability / HSA,
        HSA_score=HSA - viability,
        Bliss_CI=viability / Bliss,
        Bliss_score=Bliss - viability,
        ZIP_CI=viability / ZIP,
        ZIP_score=ZIP - viability,
        ZIP_score2=ZIP2 - viability,
        Loewe_CI=Loewe_CI,
        viability=viability,
        viability_1=viability_1,
        viability_2=viability_2,
        PANEL=unique(PANEL),
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ) ->
    combo_profiles2


combination_profiles <- combo_profiles2[,
    lapply(.SD, mean, na.rm=TRUE),
    by=c("treatment1id", "treatment2id", "sampleid", "PANEL")
]

top_synergy <- combination_profiles[
    order(-Bliss_score),
    .(
        treatment1id=unique(treatment1id)[1:10],
        treatment2id=unique(treatment2id)[1:10]
    )]

top_antagonism <- combination_profiles[
    order(Bliss_score),
    .(
        treatment1id=unique(treatment1id)[1:10],
        treatment2id=unique(treatment2id)[1:10]
    )]

subset_drugs <- rbind(top_synergy, top_antagonism)
nci_raw <- as(treatmentResponse(.nci), "data.table")
setkeyv(nci_raw, c("treatment1id", "treatment2id"))

tx1 <- unique(unname(unlist(subset_drugs)))
tx2 <- c(unique(unname(subset_drugs$treatment2id)), "")

tre <- treatmentResponse(.nci)
tre |>
    subset(treatment1id %in% tx1 & treatment2id %in% tx2, ) ->
    sub_tre

nci_dem_data <- as(sub_tre, "data.table")

## Sanity checks
combo_synergy <- treatmentResponse(nci)$profiles |>
    aggregate2(
        mean(SCORE, na.rm=TRUE),
        by=c("treatment1id", "treatment2id", "sampleid")
    )

cor(combination_profiles[,
    .(HSA_score, Bliss_score, ZIP_score, ZIP_score2)
    ], method="spearman")

cor(combination_profiles[,
    .(HSA_CI, Bliss_CI, ZIP_CI, Loewe_CI)
    ], method="spearman")

combo_plot <- melt(combination_profiles,
    measure=patterns(".*score"), value.name="score", variable.name="model")

library(ggplot2)
ggplot(combo_plot) +
    geom_density(aes(x=score, colour=model))
