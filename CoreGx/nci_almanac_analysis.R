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
        target="mono_viability",
        subset = treatment2id == "",
        viability = mean(viability) / 100,
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) |>
    endoaggregate(
        assay="mono_viability",
        target="mono_profiles",
        {
            fit <- PharmacoGx::logLogisticRegression(treatment1dose,
                                                     viability,
                                                     trunc = FALSE,
                                                     viability_as_pct = FALSE)
            ic50 <- PharmacoGx::computeIC50(treatment1dose, Hill_fit=fit)
            auc <- PharmacoGx::computeAUC(treatment1dose, Hill_fit=fit,
                area.type="Fitted")
            list(
                HS = fit[['HS']],
                E_inf = fit[['E_inf']],
                EC50 = fit[['EC50']],
                Rsquare = attributes(fit)$Rsquare,
                viability = mean(viability, na.rm = TRUE),
                auc=auc,
                ic50=ic50
            )
        },
        by=c("treatment1id", "sampleid"),
        enlist=FALSE,
        nthread=nthread
    ) -> tre
}) -> m1

# -- Attach the Hill parameters from the monotherapy profiles to the combination data
bench::system_time({
tre |>
    endoaggregate(
        assay = "sensitivity",
        target = "combo_viability",
        subset = (!is.na(treatment2dose)),
        combo_viability = mean(viability) / 100,
        by=c("treatment1id", "treatment1dose", "treatment2id", "treatment2dose",
            "sampleid")
    ) -> tre
}) -> m2

bench::system_time({
tre |>
    mergeAssays(
        "combo_viability",
        "mono_profiles",
        by = c("treatment1id", "sampleid")
    ) |>
    mergeAssays(
        "combo_viability",
        "mono_profiles",
        by.x=c("treatment2id", "sampleid"),
        by.y=c("treatment1id", "sampleid"),
        suffixes=c("_1", "_2")
    ) ->
    tre
}) -> m3
## == Vector Based computeZIPdelta ============================================
combo_profiles <- CoreGx::buildComboProfiles(ntre3, c("HS", "EC50", "E_inf", "ZIP", "combo_viability"))
bench::system_time({
combo_profiles[, 
        .computeZIPdelta(
            treatment1id = treatment1id,
            treatment2id = treatment2id,
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            sampleid = sampleid,
            HS_1 = HS_1, HS_2 = HS_2,
            EC50_1 = EC50_1, EC50_2 = EC50_2,
            E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
            combo_viability = combo_viability,
            ZIP = ZIP,
            nthread = 4,
            show_Rsqr = TRUE
        )
    ] -> delta_scores
})

# -- compute our drug synergy metrics
bench::system_time({
tre |>
    endoaggregate(
        assay="combo_viability",
        HSA = computeHSA(viability_1, viability_2),
        Bliss = computeBliss(viability_1, viability_2),
        ZIP = computeZIP(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1,
            E_inf_2 = E_inf_2,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        ),
        Loewe = PharmacoGx::computeLoewe(
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose,
            HS_1 = HS_1,
            HS_2 = HS_2,
            E_inf_1 = E_inf_1,
            E_inf_2 = E_inf_2,
            EC50_1 = EC50_1,
            EC50_2 = EC50_2
        ),
        by = assayKeys(tre, "combo_viability")
    ) ->
    ntre1
}) -> m4

bench::system_time({
ntre1 |>
    endoaggregate(
        assay = "combo_viability",
        target = "combo_scores",
        HSA_CI = combo_viability / HSA,
        HSA_score = combo_viability - HSA,
        Bliss_CI = combo_viability / Bliss,
        Bliss_score = combo_viability - Bliss,
        ZIP_CI = combo_viability / ZIP,
        ZIP_score = combo_viability - ZIP,
        Loewe_score = combo_viability - Loewe,
        Loewe_CI = loeweCI(viability = combo_viability,
                           treatment1dose = treatment1dose,
                           treatment2dose = treatment2dose,
                           HS_1 = HS_1,
                           E_inf_1 = E_inf_1,
                           EC50_1 = EC50_1,
                           HS_2 = HS_2,
                           E_inf_2 = E_inf_2,
                           EC50_2 = EC50_2),
        by = assayKeys(tre, "combo_viability"),
        nthread = nthread
    ) -> ntre2
}) -> m5

bench::system_time({ntre3 <- computeZIPdelta(
    object = ntre2,
    residual = "logcosh",
    nthread = nthread,
    show_Rsqr = TRUE
)}) -> m6
# == Analysis Ends ================================================

saveRDS(ntre3, "../../data/ntre3.rds")

drop_cols <- c("Rsquare_1", "Rsquare_2",
               "auc_1", "auc_2",
               "ic50_1", "ic50_2",
               "")
combo_viability <- tre$combo_viability

combo_viability[, (drop_cols) := NULL]

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
bench::system_time({
    combo_profiles[,
        c("viability_1", "viability_2") := .(
            PharmacoGx:::.Hill(
                log10(treatment1dose),
                c(HS_1, E_inf_1, log10(EC50_1))
            ),
            PharmacoGx:::.Hill(
                log10(treatment2dose),
                c(HS_2, E_inf_2, log10(EC50_2))
            )
        ),
        by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
    ]
})

# -- compute our drug synergy metrics
bench::system_time({
    combo_profiles |>
        aggregate2(
            HSA = PharmacoGx::computeHSA(viability_1 = viability_1,
                                       viability_2 = viability_2),
            Bliss = PharmacoGx::computeBliss(viability_1 = viability_1,
                                           viability_2 = viability_2),
            ZIP = PharmacoGx::computeZIP(
                treatment1dose = treatment1dose,
                treatment2dose = treatment2dose,
                HS_1 = HS_1,
                HS_2 = HS_2,
                EC50_1 = EC50_1,
                EC50_2 = EC50_2
            ),
            Loewe = PharmacoGx::computeLoewe(
                treatment1dose = treatment1dose,
                HS_1 = HS_1,
                E_inf_1 = E_inf_1,
                EC50_1 = EC50_1,
                treatment2dose = treatment2dose,
                HS_2 = HS_2,
                E_inf_2 = E_inf_2,
                EC50_2 = EC50_2,
                tol = 0.1
            ),
            #viability_1 = viability_1,
            #viability_2 = viability_2,
            #viability = viability / 100,
            #HS_1 = HS_1,
            #E_inf_1 = E_inf_1,
            #EC50_1 = EC50_1,
            #HS_2 = HS_2,
            #E_inf_2 = E_inf_2,
            #EC50_2 = EC50_2,
            nthread = 1,
            by = c("treatment1id", "treatment2id",
                   "treatment1dose","treatment2dose", "sampleid")
        ) -> combo_viability
})
treatmentResponse(nci)$combo_viability <- combo_viability
profiles <- c("HS", "E_inf", "EC50")
j <- parse(text = paste0(".(",
        paste(c("treatment1id", "sampleid", profiles), collapse = ","),
    ")"))
combo_keys <- c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
monotherapy_params <- treatmentResponse(nci)$monotherapy_profiles[,
    eval(j)
]
j <- parse(text = paste0(".(",
        paste(c(combo_keys, profiles), collapse = ","),
    ")"))
combo_dose <- combo_dose[
    treatmentResponse(nci)$monotherapy_profiles,
    eval(j),
    on = c(treatment1id = "treatment1id",
           sampleid = "sampleid")
]
combo_profiles_cols <- c(
    paste0(profiles, sep = "_1"),
    paste0(profiles, sep = "_2")
)
j <- parse(text = paste0(".(",
        paste(c(combo_keys, combo_profiles_cols), collapse = ","),
    ")"))
combo_dose <- merge(
    combo_dose,
    monotherapy_params,
    by.x=c("treatment2id", "sampleid"),
    by.y=c("treatment1id", "sampleid"),
    suffixes=c("_1", "_2")
)[, eval(j)]
setkeyv(combo_dose, c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid"))
delta_score_input <- combo_dose[treatmentResponse(nci)$combo_viability][,
                                    c(colnames(combo_dose), "ZIP", "viability"),
                                    with = FALSE
                                ]

delta_score_input[delta_score_input[.N, treatment1id]]

bench::system_time({
fit_1_to_2 <- delta_score_input |>
    aggregate(
        treatment1id = treatment1id,
        treatment2id = treatment2id,
        treatment1dose = treatment1dose,
        treatment2dose = treatment2dose,
        sampleid = sampleid,
        EC50_1 = EC50_1,
        HS_1 = HS_1,
        viability = viability,
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid")
    ) |>
    aggregate(
        PharmacoGx::estimateNewPotency(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1)
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    )
})

bench::system_time({
fit_1_to_2 <- delta_score_input[treatment1id == delta_score_input[.N, treatment1id] & treatment2id == delta_score_input[.N, treatment2id] &
#sampleid == "KM12"
viability < 1
] |>
    aggregate(
        PharmacoGx::estimateNewPotency(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1)
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE,
        nthread = 3
    )
})


bench::system_time({
    combo_viability |>
        aggregate2(
            Loewe = Loewe,
            Loewe_CI1 = PharmacoGx::LoeweCI(viability = Loewe,
                                            treatment1dose = treatment1dose,
                                            treatment2dose = treatment2dose,
                                            HS_1 = HS_1, HS_2 = HS_2,
                                            E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
                                            EC50_1 = EC50_1, EC50_2 = EC50_2),
            viability_1 = viability_1,
            viability_2 = viability_2,
            viability = viability / 100,
            HS_1 = HS_1,
            E_inf_1 = E_inf_1,
            EC50_1 = EC50_1,
            HS_2 = HS_2,
            E_inf_2 = E_inf_2,
            EC50_2 = EC50_2,
            nthread = 1,
            by = c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid")
        ) -> Loewe_check
})
tol <- 1e-1
Loewe_check[abs(Loewe_CI1 - 1) > tol, hist(x = Loewe_CI1)]

object <- treatmentResponse(nci)
combo_viability <- object$combo_viability

bench::system_time({
    combo_profiles <- buildComboProfiles(object, c("HS", "EC50", "E_inf", "ZIP", "viability"))
})

# -- compute combination index and score
bench::system_time({
    combo_viability |>
        aggregate2(
            HSA_CI = viability / HSA,
            HSA_score = HSA - viability,
            Bliss_CI = viability / Bliss,
            Bliss_score = Bliss - viability,
            ZIP_CI = viability / ZIP,
            ZIP_score = ZIP - viability,
            Loewe_score = Loewe - viability,
            by=c("treatment1id", "treatment2id",
                 "treatment1dose", "treatment2dose",
                 "sampleid"),
            nthread = 1
        ) -> combo_scores
})

dt <- merge.data.table(combo_profiles, combo_profiles1, by=c("treatment1id", "treatment2id", "treatment1dose", "treatment2dose", "sampleid"))

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

test_combo <- combo_profiles[treatment1id == "Cyclophosphamide" & treatment2id == "Vorinostat" & sampleid == "RPMI-8226"]

result <- PharmacoGx::projectedResponse(
    dose_add = test_combo[, treatment1dose],
    EC50_add = test_combo[, EC50_1],
    HS_add = test_combo[, HS_1],
    dose_to = test_combo[, treatment2dose],
    EC50_proj = res[["EC50"]],
    HS_proj = res[["HS"]]
)


treatmentResponse(nci) |>
    subset(treatment2id == "") |>
    aggregate(
        assay="sensitivity",
        viability=mean(viability),
        by=c("treatment1id", "treatment1dose", "sampleid")
    ) -> mono_viability

mono_profiles <- monotherapy_profiles[!is.na(Rsq)]
temp_monotherapy_profiles <- mono_viability[
    mono_profiles,,
    on = c(treatment1id = "treatment1id",
           sampleid = "sampleid")
]

tre |>
    subset(treatment1id %in% select_treatment1id &
           treatment2id %in% select_treatment2id,
           sampleid %in% select_sampleid
    ) -> ntre_small
## test computeZIPdelta ======================
object <- treatmentResponse(nci)
bench::system_time({
    combo_profiles <- buildComboProfiles(object, c("HS", "EC50", "E_inf", "ZIP", "viability"))
})
combo_profiles[, E_inf_1 := E_inf_1 / 100]
combo_profiles[, E_inf_2 := E_inf_2 / 100]

combo_viability <- object$combo_viability
combo_scores <- object$combo_scores
## subset data ===============================
select_treatment1id <- c("1,3-Bis(2-chloroethyl)-1-nitrosourea", "Zolendronic Acid")
select_treatment2id <- c("2-fluoroAraA (fludarabine)", "Vorinostat")
select_sampleid <- c("786-O", "UO-31")
combo_profiles <- combo_profiles[
    treatment1id %in% select_treatment1id &
    treatment2id %in% select_treatment2id &
    sampleid %in% select_sampleid,
]
combo_viability <- combo_viability[
    treatment1id %in% select_treatment1id &
    treatment2id %in% select_treatment2id &
    sampleid %in% select_sampleid,
]
combo_scores <- combo_scores[
    treatment1id %in% select_treatment1id &
    treatment2id %in% select_treatment2id &
    sampleid %in% select_sampleid,
]
setkeyv(combo_profiles, combo_keys)
# args =======================================

bench::system_time({
    combo_twowayFit <- fitTwowayZIP(combo_profiles, residual = "logcosh",
                                    show_Rsqr = TRUE, nthread = nthread)
})

combo_profiles |>
    aggregate(
        estimateProjParams(
            dose_to = treatment1dose,
            viability = viability,
            dose_add = unique(treatment2dose),
            EC50_add = unique(EC50_2),
            HS_add = unique(HS_2),
            use_L2 = FALSE,
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        enlist = FALSE
    ) -> fit_2_to_1
## examine varying dose vs. viability 
combo_profiles[, .(treatment1id, treatment2id, sampleid, treatment1dose, viability),
               by = .(treatment1id, treatment2id, sampleid, treatment2dose)]
combo_profiles |>
    aggregate(
        estimateNewPotency(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1),
            use_L2 = FALSE,
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    ) -> fit_1_to_2
# L2 ====================
combo_profiles |>
    aggregate(
        estimateNewPotency(
            dose_to = treatment1dose,
            viability = viability,
            dose_add = unique(treatment2dose),
            EC50_add = unique(EC50_2),
            HS_add = unique(HS_2),
            use_L2 = TRUE,
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        enlist = FALSE
    ) -> fit_2_to_1_L2
combo_profiles |>
    aggregate(
        estimateNewPotency(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1),
            use_L2 = TRUE,
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    ) -> fit_1_to_2_L2

combo_profiles <- combo_profiles[
    fit_1_to_2, ,
    on = c(
        treatment1id = "treatment1id",
        treatment2id = "treatment2id",
        treatment1dose = "treatment1dose",
        sampleid = "sampleid"
    )
]


combo_profiles <- merge(
    combo_profiles,
    fit_2_to_1,
    by.x=c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    by.y=c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    suffixes=c("_1_to_2", "_2_to_1")
)

combo_profiles |>
    aggregate(
        delta_score = .computeZIPDelta(
            EC50_1_to_2 = EC50_proj_1_to_2,
            EC50_2_to_1 = EC50_proj_2_to_1,
            EC50_1 = EC50_1, EC50_2 = EC50_2,
            HS_1_to_2 = HS_proj_1_to_2,
            HS_2_to_1 = HS_proj_2_to_1,
            HS_1 = HS_1, HS_2 = HS_2,
            E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
            treatment1dose = treatment1dose,
            treatment2dose = treatment2dose
        ),
        by = combo_keys,
        nthread = nthread
    ) -> delta_scores
# ============================================
bench::system_time({
    combo_viability_small |>
        aggregate2(
            HSA_CI = viability / HSA,
            HSA_score = HSA - viability,
            Bliss_CI = viability / Bliss,
            Bliss_score = Bliss - viability,
            ZIP_CI = viability / ZIP,
            ZIP_score = ZIP - viability,
            Loewe_score = Loewe - viability,
            by=c("treatment1id", "treatment2id",
                 "treatment1dose", "treatment2dose",
                 "sampleid"),
            nthread = 1
        ) -> combo_scores
})
combo_scores <- combo_scores[delta_scores,,
    on = c(treatment1id = "treatment1id",
           treatment2id = "treatment2id",
           treatment1dose = "treatment1dose",
           treatment2dose = "treatment2dose",
           sampleid = "sampleid")]

res <- PharmacoGx::estimateNewPotency(
    viability = test_combo[1:3, viability],
    dose_add = test_combo[1, treatment1dose],
    EC50_add = test_combo[1, EC50_1],
    HS_add = test_combo[1, HS_1],
    dose_to = test_combo[1:3, treatment2dose]
)
res_L2 <- PharmacoGx::estimateNewPotency(
    viability = test_combo[1:3, viability],
    dose_add = test_combo[1, treatment1dose],
    EC50_add = test_combo[1, EC50_1],
    HS_add = test_combo[1, HS_1],
    dose_to = test_combo[1:3, treatment2dose],
    use_L2 = TRUE
)

test_combo <- combo_profiles[treatment1id == "Cyclophosphamide" & treatment2id == "Vorinostat" & sampleid == "RPMI-8226"]

test_combo_L2 <- copy(test_combo)

bench::system_time({
test_combo_normal |>
    aggregate(
        estimateProjParams(
            dose_to = treatment1dose,
            viability = viability,
            dose_add = unique(treatment2dose),
            EC50_add = unique(EC50_2),
            HS_add = unique(HS_2),
            E_inf_add = unique(E_inf_2),
            residual = residual,
            show_Rsqr = show_Rsqr
        ),
        moreArgs = list( 
            residual = "normal",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        enlist = FALSE
    ) -> fit_2_to_1_test_normal

test_combo_normal |>
    aggregate(
        estimateProjParams(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1),
            E_inf_add = unique(E_inf_1),
            residual = residual,
            show_Rsqr = show_Rsqr
        ),
        moreArgs = list(
            residual = "normal",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    ) -> fit_1_to_2_test_normal

test_combo_normal <- test_combo_normal[
    fit_1_to_2_test_normal, ,
    on = c(
        treatment1id = "treatment1id",
        treatment2id = "treatment2id",
        treatment1dose = "treatment1dose",
        sampleid = "sampleid"
    )
]
    
test_combo_normal <- merge(
    test_combo_normal,
    fit_2_to_1_test_normal,
    by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    suffixes = c("_1_to_2", "_2_to_1")
)
})

bench::system_time({
test_combo_cauchy |>
    aggregate(
        estimateProjParams(
            dose_to = treatment1dose,
            viability = viability,
            dose_add = unique(treatment2dose),
            EC50_add = unique(EC50_2),
            HS_add = unique(HS_2),
            E_inf_add = unique(E_inf_2),
            residual = residual,
            show_Rsqr = show_Rsqr
        ),
        moreArgs = list( 
            residual = "Cauchy",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        enlist = FALSE
    ) -> fit_2_to_1_test_cauchy

test_combo_cauchy |>
    aggregate(
        estimateProjParams(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1),
            E_inf_add = unique(E_inf_1),
            residual = residual,
            show_Rsqr = show_Rsqr
        ),
        moreArgs = list(
            residual = "Cauchy",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    ) -> fit_1_to_2_test_cauchy

test_combo_cauchy <- test_combo_cauchy[
    fit_1_to_2_test_cauchy, ,
    on = c(
        treatment1id = "treatment1id",
        treatment2id = "treatment2id",
        treatment1dose = "treatment1dose",
        sampleid = "sampleid"
    )
]
    
test_combo_cauchy <- merge(
    test_combo_cauchy,
    fit_2_to_1_test_cauchy,
    by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    suffixes = c("_1_to_2", "_2_to_1")
)
})

# =======
test_combo_default <- copy(test_combo)

test_combo <- test_combo_L2

test_combo <- test_combo_drc

unique_t1_dose <- unique(test_combo[, treatment1dose])
unique_t2_dose <- unique(test_combo[, treatment2dose])
cols <- palette(rainbow(3))
plot(NULL, xlim = c(-10, 10), ylim = c(0, 2),
     ylab = "Viability of adding Vorinostat to Cyclophosphamide (RPMI-8226)",
     xlab = "log([Cyclophosphamide])")
for (i in seq_along(unique_t2_dose)) {
    dose_add <- unique_t2_dose[i]
    EC50_proj <- unique(test_combo[treatment2dose == dose_add, EC50_proj_2_to_1])
    HS_proj <- unique(test_combo[treatment2dose == dose_add, HS_proj_2_to_1])
    EC50_add <- unique(test_combo[treatment2dose == dose_add, EC50_2])
    HS_add <- unique(test_combo[treatment2dose == dose_add, HS_2])
    E_inf_add <- unique(test_combo[treatment2dose == dose_add, E_inf_2])
    E_inf_proj <- unique(test_combo[treatment2dose == dose_add, E_inf_proj_2_to_1])
    dose_to <- test_combo[treatment2dose == dose_add, treatment1dose]
    y <- test_combo[treatment2dose == dose_add, viability]
    E_min_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
    curve(Hill_4par(
            E_nnf = E_min_proj,
            E_inf = E_inf_proj,
            #E_inf = 0,
            HS = HS_proj,
            EC50 = log10(EC50_proj),
            dose = x
            ),
          from = -6, to = 6, add = TRUE,
          col = cols[i])
    points(x = log10(dose_to), y = y, col = cols[i])
}
legend(-6, 0.75, 
    legend = c("[Vorinostat] = 0.01", "[Vorinostat] = 0.10", "[Vorinostat] = 1.00"), 
    col = cols, 
    lty = 1,
    box.lty = 0
)

opar <- par()

vorinostat_profiles <- object$monotherapy_profiles[treatment1id == "Vorinostat" & sampleid == "RPMI-8226"]

plotProjHill <- function(combo_twowayFit, treatment1, treatment2, cellline, title = NULL) {
    select_combo <- combo_twowayFit[treatment1id == treatment1 &
                                    treatment2id == treatment2 &
                                    sampleid == cellline]
    if (dim(select_combo)[1] <= 0)
        stop("no combo")
    unique_t1_dose <- unique(select_combo[, treatment1dose])
    unique_t2_dose <- unique(select_combo[, treatment2dose])
    Rsqr_1_to_2 <- vector(mode = "numeric", length = length(unique_t1_dose))
    Rsqr_2_to_1 <- vector(mode = "numeric", length = length(unique_t2_dose))

    cols <- palette(rainbow(length(unique_t1_dose)))

    if (is.null(title))
        title <- deparse(substitute(combo_twowayFit))


    plot_2_to_1 <- ggplot(data = select_combo) +
        xlim(-10, 10) +
        ylim(0, 2 * max(select_combo$viability)) +
        labs(
            x = paste0("log10([", treatment2, "])"), y = paste(cellline, "viability"),
            color = paste0("[", treatment1, "]")
        ) +
        geom_point(aes(
            x = log10(treatment1dose),
            y = viability,
            colour = as.factor(treatment2dose)
        ))
    for (i in seq_along(unique_t2_dose)) {
        dose_add <- unique_t2_dose[i]
        EC50_proj <- unique(select_combo[treatment2dose == dose_add, EC50_proj_2_to_1])
        HS_proj <- unique(select_combo[treatment2dose == dose_add, HS_proj_2_to_1])
        EC50_add <- unique(select_combo[treatment2dose == dose_add, EC50_2])
        E_inf_add <- unique(select_combo[treatment2dose == dose_add, E_inf_2])
        E_inf_proj <- unique(select_combo[treatment2dose == dose_add, E_inf_proj_2_to_1])
        HS_add <- unique(select_combo[treatment2dose == dose_add, HS_2])
        dose_to <- select_combo[treatment2dose == dose_add, treatment1dose]
        E_nnf_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
        Rsqr_2_to_1[i] <- unique(select_combo[treatment2dose == dose_add, Rsqr_2_to_1])
        plot_2_to_1 <- plot_2_to_1 + geom_function(fun = Hill_4par,
            args = list(
                HS = HS_proj,
                E_nnf = E_nnf_proj,
                E_inf = E_inf_proj,
                EC50 = EC50_proj
            ),
            colour = cols[i]
        )
    }
    plot_2_to_1 <- plot_2_to_1 +
        scale_colour_manual(values = cols) +
        theme(legend.position = "bottom") +
        annotate(geom = "text",
                 label = paste("R^2", "=" , round(Rsqr_2_to_1, 4)),
                 x = seq(-6, 6, length.out = length(Rsqr_2_to_1)),
                 y = 1.5, colour = cols, size = 4)

    legend(-10, 0.5, 
        legend = paste0("[", treatment2, "] = ", unique_t2_dose,
                        ", R square = ", round(Rsqr_2_to_1, digits = 4)), 
        col = cols, 
        lty = 1,
        box.lty = 0
    )

    plot(
        NULL, xlim = c(-10, 10), ylim = c(0, 2),
        ylab = paste("Response of adding", treatment1, "to", treatment2),
        xlab = paste0("log10([", treatment2,"])"),
        main = title
    )
    for (i in seq_along(unique_t1_dose)) {
        dose_add <- unique_t1_dose[i]
        EC50_proj <- unique(select_combo[treatment1dose == dose_add, EC50_proj_1_to_2])
        HS_proj <- unique(select_combo[treatment1dose == dose_add, HS_proj_1_to_2])
        EC50_add <- unique(select_combo[treatment1dose == dose_add, EC50_1])
        HS_add <- unique(select_combo[treatment1dose == dose_add, HS_1])
        E_inf_add <- unique(select_combo[treatment1dose == dose_add, E_inf_1])
        E_nnf_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
        E_inf_proj <- unique(select_combo[treatment1dose == dose_add, E_inf_proj_1_to_2])
        dose_to <- select_combo[treatment1dose == dose_add, treatment2dose]
        Rsqr_1_to_2[i] <- unique(select_combo[treatment1dose == dose_add, Rsqr_1_to_2])
        y <- select_combo[treatment1dose == dose_add, viability]
        curve(
            PharmacoGx::Hill_4par(
                E_nnf = E_nnf_proj,
                E_inf = E_inf_proj,
                HS = HS_proj,
                EC50 = log10(EC50_proj),
                dose = x
            ),
            from = -10, to = 10, add = TRUE, col = cols[i]
        )
        points(x = log10(dose_to), y = y, col = cols[i])
    }
    legend(-10, 0.5, 
        legend = paste0("[", treatment1, "] = ", unique_t1_dose,
                        ", R square = ", round(Rsqr_1_to_2, digits = 4)), 
        col = cols, 
        lty = 1,
        box.lty = 0
    )
}

.plotProjHill(test_combo_logcosh2,
             treatment1 = "Methotrexate",
             treatment2 = "Zolendronic Acid",
             cellline = "UO-31")

combo_viability <- ntre3$combo_viability

demo_combo <- combo_viability[
    treatment1id == "Idarubicin" &
    treatment2id == "Erlotinib" &
    sampleid == "OVCAR-4"
]

select_combo <- fitTwowayZIP(demo_combo)
treatment1 <- "Idarubicin"
treatment2 <- "Erlotinib"
plot(
    NULL, xlim = c(-6, 6), ylim = c(0, 1),
    ylab = paste("Viability of OVCAR-4"),
    xlab = paste0("log10([", treatment2,"])"),
    main = "Dose-effect Curve of Adding Idarubicin to Erlotinib"
)
dose_add <- 1e-02
EC50_proj <- unique(select_combo[treatment1dose == dose_add, EC50_proj_1_to_2])
HS_proj <- unique(select_combo[treatment1dose == dose_add, HS_proj_1_to_2])
EC50_add <- unique(select_combo[treatment1dose == dose_add, EC50_1])
E_inf_add <- unique(select_combo[treatment1dose == dose_add, E_inf_1])
E_inf_proj <- unique(select_combo[treatment1dose == dose_add, E_inf_proj_1_to_2])
HS_add <- unique(select_combo[treatment1dose == dose_add, HS_1])
dose_to <- select_combo[treatment1dose == dose_add, treatment2dose]
E_ninf_proj <- PharmacoGx:::.Hill(log10(dose_add), c(HS_add, E_inf_add, log10(EC50_add)))
Rsqr_1_to_2 <- unique(select_combo[treatment1dose == dose_add, Rsqr_1_to_2])
y <- select_combo[treatment1dose == dose_add, combo_viability]
curve(
    PharmacoGx::hill4Par(
        E_ninf = E_ninf_proj,
        E_inf = E_inf_proj,
        HS = HS_proj,
        EC50 = log10(EC50_proj),
        dose = x
    ),
    from = -10, to = 10, add = TRUE
    )
points(x = log10(dose_to), y = y)
legend(-6.5, 0.6,
    legend = paste0("[", treatment2, "] = ", dose_add,
                    ", R square = ", round(Rsqr_1_to_2, digits = 4)),
    lty = 1,
    box.lty = 0
)
abline(h = E_inf_proj, col = "red", lty = 2)
text(5, E_inf_proj - 0.05, labels = expression("E"[infinity]), col = "red", cex = 1.5)
abline(h = E_ninf_proj, col = "blue", lty = 2)
text(-5, E_ninf_proj - 0.05, labels = expression("E"[-infinity]), col = "blue", cex = 1.5)
abline(v = log10(EC50_proj), col = "green", lty = 2)
text(log10(EC50_proj), 0, labels = "EC50", col = "green", cex = 1.5)

## examine varying dose vs. viability 
test_combo[, .(treatment1id, treatment2id, sampleid, treatment1dose, viability),
           by = .(treatment1id, treatment2id, sampleid, treatment2dose)]

bench::system_time({test_combo_logcosh_twowayFit <- fitTwowayZIP(
    test_combo_logcosh,
    residual = "normal",
    show_Rsqr = TRUE)})

bench::system_time({
 test_combo_logcosh |>
     aggregate(
         delta_score = .computeZIPdelta(
             treatment1id = treatment1id,
             treatment2id = treatment2id,
             treatment1dose = treatment1dose,
             treatment2dose = treatment2dose,
             sampleid = sampleid,
             viability = viability,
             HS_1 = HS_1, HS_2 = HS_2,
             EC50_1 = EC50_1, EC50_2 = EC50_2,
             E_inf_1 = E_inf_1, E_inf_2 = E_inf_2,
             residual = residual
         ),
         treatment1dose = treatment1dose,
         treatment2dose = treatment2dose,
         moreArgs = list(residual = "logcosh"),
         by = c("treatment1id", "treatment2id", "sampleid")
     ) -> test_combo_logcosh_delta
})

test_combo[,.(Rsqr_2_to_1, Rsqr_1_to_2),
           by = .(treatment1id, treatment2id, sampleid, treatment1dose, treatment2dose)]
test_combo_L2[,.(Rsqr_2_to_1, Rsqr_1_to_2),
              by = .(treatment1id, treatment2id, sampleid, treatment1dose, treatment2dose)]

# ========================================================================
bench::system_time({
    combo_profiles <- buildComboProfiles(ntre3, c("HS", "EC50", "E_inf", "ZIP", "combo_viability"))
})
hard_combo <- combo_profiles[treatment1id == "Zolendronic Acid" &
                             treatment2id == "2-fluoroAraA (fludarabine)" &
                             sampleid == "786-O"]

hard_combo <- combo_profiles[treatment2id == "Zolendronic Acid" &
                             treatment1id == "Methotrexate" &
                             sampleid == "UO-31"]

combo_fit2Way <- fitTwowayZIP(hard_combo, show_Rsqr = TRUE)

.plotProjHill(combo_fit2Way)
# ========================================================================


test_combo_drc <- copy(test_combo)

bench::system_time({
test_combo_drc |>
    aggregate(
        estimateProjParams(
            dose_to = treatment1dose,
            viability = viability,
            dose_add = unique(treatment2dose),
            EC50_add = unique(EC50_2),
            HS_add = unique(HS_2),
            E_inf_add = unique(E_inf_2),
            show_Rsqr = show_Rsqr,
            residual = residual
        ),
        moreArgs = list(
            residual = "drc",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
        enlist = FALSE
    ) -> fit_2_to_1_test_drc
test_combo_drc |>
    aggregate(
        estimateProjParams(
            dose_to = treatment2dose,
            viability = viability,
            dose_add = unique(treatment1dose),
            EC50_add = unique(EC50_1),
            HS_add = unique(HS_1),
            E_inf_add = unique(E_inf_1),
            residual = residual,
            show_Rsqr = TRUE
        ),
        moreArgs = list(
            residual = "drc",
            show_Rsqr = TRUE
        ),
        by = c("treatment1id", "treatment2id", "treatment1dose", "sampleid"),
        enlist = FALSE
    ) -> fit_1_to_2_test_drc

test_combo_drc <- test_combo_drc[
    fit_1_to_2_test_drc, ,
    on = c(
        treatment1id = "treatment1id",
        treatment2id = "treatment2id",
        treatment1dose = "treatment1dose",
        sampleid = "sampleid"
    )
]
    
test_combo_drc <- merge(
    test_combo_drc,
    fit_2_to_1_test_drc,
    by.x = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    by.y = c("treatment1id", "treatment2id", "treatment2dose", "sampleid"),
    suffixes = c("_1_to_2", "_2_to_1")
)
})


test_combo[,.(Rsqr_2_to_1, Rsqr_1_to_2),
           by = .(treatment1id, treatment2id, sampleid, treatment1dose, treatment2dose)]
test_combo_L2[,.(Rsqr_2_to_1, Rsqr_1_to_2),
              by = .(treatment1id, treatment2id, sampleid, treatment1dose, treatment2dose)]

# == computeHSA ==========================================================


