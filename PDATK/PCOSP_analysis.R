library(MetaGxPancreas)
library(PDATK)
library(BiocParallel)
library(msigdbr)
library(ggplot2)
library(data.table)
library(qs)


## ---- 0. Configure script parameters


## WARNING: This script requires > 32 GB of RAM to run
randomSeed <- 1987
numModels <- 1000
minDaysSurvived <- 365


## ---- 1. Pancreatic Cancer Overall Survival Predictor


## -- 1.0 Download and preprocess the data

# NOTE: This will download several GiB of data to your 
#>ExperimentHub Cache
metaGxPanc <- loadPancreasDatasets()
pancreasCohorts <- metaGxPanc$SummarizedExperiments

cohortsWithSurvival <- c("ICGCMICRO", "ICGCSEQ", "PCSI", 
    "TCGA", "KIRBY", "UNC", "CHEN", "COLLISON", "ZHANG", 
    "OUH", "WINTER")
names(pancreasCohorts) <- gsub('_SumExp', '', names(pancreasCohorts))
pancreasCohortsWithSurv <- pancreasCohorts[cohortsWithSurvival]

# SurvivalExperiment is just a wrapper around SummarizedExpeirment
pancreasSurvExps <- lapply(pancreasCohortsWithSurv, FUN=SurvivalExperiment, 
    survival_time='days_to_death', event_occurred='vital_status')
# CohortList is a SimpleList of SurvivalExperiments
pancreasCohortList <- CohortList(pancreasSurvExps)

## -- 1.1 Split Training and Validation Data

# Find common genes
commonGenes <- findCommonGenes(pancreasCohortList)

# subset on a CohortList subsets all items
cohortList <- subset(pancreasCohortList, subset=commonGenes)

# Train from the ICGC data only
ICGCcohorts <- grepl('ICGC', names(cohortList))
ICGCcohortList <- cohortList[ICGCcohorts]
# The remainder of the data will be used for validation
validationCohortList <- cohortList[!ICGCcohorts]

validationCohortList <- dropNotCensored(validationCohortList)
ICGCcohortList <- dropNotCensored(ICGCcohortList)

# Find common samples between the sequencing and array data from ICGC
commonSamples <- findCommonSamples(ICGCcohortList)

# split into shared samples for training, the rest for testing
ICGCtrainCohorts <- subset(ICGCcohortList, select=commonSamples)
ICGCtestCohorts <- subset(ICGCcohortList, select=commonSamples, 
    invert=TRUE)

# merge our training cohort test data into the rest of the validation
#>data
validationCohortList <- c(ICGCtestCohorts, validationCohortList)

# drop ICGCSEQ from the validation data, because it only has 7 
#>patients
validationCohortList <- 
    validationCohortList[names(validationCohortList) != 'ICGCSEQ']

## -- 1.2 Setup a PCOSP model

# Set a seed for a reproducible analysis!
set.seed(randomSeed)

# Construct a PCOSP model object
PCOSPmodel <- PCOSP(ICGCtrainCohorts, minDaySurvived=minDaysSurvived,
    randomSeed=randomSeed)

# -- 1.3 Training a PCOSP model

trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=numModels, 
    minAccuracy=0.6)

qsave(trainedPCOSPmodel, file='results/trainedPCOSP.qs')

# -- 1.4 Risk Prediction with a PCOSP model

PCOSPpredValCohorts <- predictClasses(validationCohortList,
    model=trainedPCOSPmodel)

mcols(PCOSPpredValCohorts)

## -- 1.5 Validating a PCOSP model

validatedPCOSPmodel <- validateModel(trainedPCOSPmodel,
    valData=PCOSPpredValCohorts)

qsave(validatedPCOSPmodel, file='results/validatedPCOSP.qs')

## -- 1.6 Plotting Model Performance

PCOSPdIndexForestPlot <- forestPlot(validatedPCOSPmodel, 
    stat='log_D_index')

ggsave(PCOSPdIndexForestPlot, width=8.5, height=11,
    file='results/figures/PCOSPlogDindexForestPlot.pdf')

PCOSPconcIndexForestPlot <- forestPlot(validatedPCOSPmodel, 
    stat='concordance_index')

ggsave(PCOSPconcIndexForestPlot, width=8.5, height=11,
    file='results/figures/PCOSPcindexForestPlot.pdf')

cohortROCplots <- plotROC(validatedPCOSPmodel, alpha=0.05)

ggsave(cohortROCplots, width=8, height=8,
    file='results/figures/PCOSPcohortROC.pdf')


## ---- 2. Permutation Testing


## -- 2.1 Random Label Shuffling

# Merge to a single SurvivalExperiment
ICGCtrainCohorts <- merge(ICGCtrainCohorts[[1]], 
    ICGCtrainCohorts[[2]], cohortNames=names(ICGCtrainCohorts))
RLSmodel <- RLSModel(ICGCtrainCohorts, minDaysSurvived=minDaysSurvived, 
    randomSeed=randomSeed)

trainedRLSmodel <- trainModel(RLSmodel, numModels=numModels)

qsave(trainedRLSmodel, file='results/trainedRLS.qs')

RLSpredCohortList <- predictClasses(validationCohortList, 
    model=trainedRLSmodel)

validatedRLSmodel <- validateModel(trainedRLSmodel, 
    RLSpredCohortList)

qsave(validatedRLSmodel, file='results/validatedRLS.qs')

# Plots
RLSmodelComparisonPlot <- densityPlotModelComparison(validatedRLSmodel,
    validatedPCOSPmodel, title='Random Label Shuffling vs PCOSP',
    mDataTypeLabels=c(rna_seq='Sequencing-based', rna_micro='Array-based',
        combined='Overall'))
ggsave(RLSmodelComparisonPlot, width=8, height=6, 
    file='results/figures/RLSmodelComparisonDensity.pdf')

## -- 2.2 Random Gene Assignment

RGAmodel <- RGAModel(ICGCtrainCohorts, randomSeed=randomSeed)

trainedRGAmodel <- trainModel(RGAmodel, numModels=numModels)

qsave(trainedRGAmodel, file='results/trainedRGA.qs')

RGApredCohortList <- predictClasses(validationCohortList,
    model=trainedRGAmodel)

validatedRGAmodel <- validateModel(trainedRGAmodel, 
    RGApredCohortList)

qsave(validatedRGAmodel, file='results/validatedRGA.qs')

# Plots
RGAmodelComparisonPlot <-  densityPlotModelComparison(validatedRGAmodel,
    validatedPCOSPmodel, title='Random Gene Assignment vs PCOSP',
    mDataTypeLabels=c(rna_seq='Sequencing-based', rna_micro='Array-based',
        combined='Overall'))

ggsave(RGAmodelComparisonPlot, width=8, height=6, 
    file='results/figures/RGAmodelComparisonDensity.pdf')


## ---- 3. Pathway Analysis


## -- 3.1 Get Top Predictive Featuress
topFeatures <- getTopFeatures(validatedPCOSPmodel)

## -- 3.2 Query Genesets for Enriched Pathways
allHumanGeneSets <- msigdbr()
allGeneSets <- as.data.table(allHumanGeneSets)
geneSets <- allGeneSets[grepl('^GO.*|.*CANONICAL.*|^HALLMARK.*', gs_name),
    .(gene_symbol, gs_name)]

GSEAresultDT <- runGSEA(validatedPCOSPmodel, geneSets)
fwrite(GSEAresultDT, file='results/GSEAresultDT.csv')


## ---- 4. Clinical Model Comparison


## -- 4.1 Build the Model
clinicalModel <- ClinicalModel(ICGCtrainCohorts,
    formula='prognosis ~ sex + age + T + N + M + grade',
    randomSeed=randomSeed)

## -- 4.2 Train the Model
trainedClinicalModel <- trainModel(clinicalModel)

qsave(trainedClinicalModel, file='results/trainedClinical.qs')

## -- 4.3 Predict the Classes
hasModelParamsCohortList <-
    PCOSPpredValCohorts[c('ICGCMICRO', 'TCGA', 'PCSI', 'OUH')]

clinicalPredCohortList <- predictClasses(hasModelParamsCohortList,
    model=trainedClinicalModel)

## -- 4.4 Validate the Model
validatedClinicalModel <- validateModel(trainedClinicalModel,
    clinicalPredCohortList)

qsave(validatedClinicalModel, file='results/validatedClinical.qs')

## -- 4.5 Visualize the Comparsion with PCOSP
clinicalVsPCOSPbarPlot <- barPlotModelComparison(validatedClinicalModel,
    validatedPCOSPmodel, stat='AUC')
ggsave(clinicalVsPCOSPbarPlot, width=8, height=8,
    file='results/figures/clinicalVsPCOSPbarPlot.pdf')

clinicalVsPCOSP <- compareModels(validatedClinicalModel, validatedPCOSPmodel)

clinVsPCOSPdIndexForestPlot <- forestPlot(clinicalVsPCOSP, stat='log_D_index')

ggsave(clinVsPCOSPdIndexForestPlot, width=8.5, height=11,
    file='results/figures/clinVsPCOSPdIndexForestPlot.pdf')

clinVsPCOSPconcIndexForestPlot <- forestPlot(clinicalVsPCOSP,
    stat='concordance_index')

ggsave(clinVsPCOSPconcIndexForestPlot, width=8.5, height=11,
    file='results/figures/clinVsPCOSPconcIndexForestPlot.pdf')


## ---- 5 Existing Classifier Comparison


## -- 5.1 Make the models
chenGeneFuModel <- GeneFuModel(randomSeed=randomSeed)
birnbaumGeneFuModel <- GeneFuModel(randomSeed=randomSeed)
haiderGeneFuModel <- GeneFuModel(randomSeed=randomSeed)

# load the classifier centroids
data(existingClassifierData)

models(chenGeneFuModel) <- SimpleList(list(chen=chen))
models(birnbaumGeneFuModel) <- SimpleList(list(birnbuam=birnbaum))
models(haiderGeneFuModel) <- SimpleList(list(haider=NA)) 

## --- 5.2 Predict the Classes

chenClassPredictions <- 
    predictClasses(PCOSPpredValCohorts[names(haiderSigScores)],
        model=chenGeneFuModel)
birnClassPredictions <- 
    predictClasses(PCOSPpredValCohorts[names(haiderSigScores)],
        model=birnbaumGeneFuModel)

haiderClassPredictions <- 
    PCOSPpredValCohorts[names(haiderSigScores)]
# Manually assign the scores to the prediction cohorts
for (i in seq_along(haiderClassPredictions)) {
    colMData <- colData(haiderClassPredictions[[i]])
    colMData$genefu_score <- NA_real_
    colMData[
        rownames(colMData) %in% names(haiderSigScores[[i]]), 
        ]$genefu_score <- 
            haiderSigScores[[i]][
                names(haiderSigScores[[i]]) %in% 
                    rownames(colMData)
            ]
    colData(haiderClassPredictions[[i]]) <- colMData
}
# Setup the correct model metadata
mcols(haiderClassPredictions)$hasPredictions <- TRUE
metadata(haiderClassPredictions)$predictionModel <- haiderGeneFuModel

## -- 5.3 Validate the Models
validatedChenModel <- validateModel(chenGeneFuModel, 
    valData=chenClassPredictions)
validatedBirnbaumModel <- validateModel(birnbaumGeneFuModel, 
    valData=birnClassPredictions)
validatedHaiderModel <- validateModel(haiderGeneFuModel, 
    valData=haiderClassPredictions)

qsave(validatedChenModel, file='results/validatedChen.qs')
qsave(validatedBirnbaumModel, file='results/validatedBirnbaum.qs')
qsave(validatedHaiderModel, file='results/validatedHaider.qs')

## -- 5.4 Model Performance Meta-Analysis

genefuModelComparisons <- compareModels(validatedChenModel,
    validatedBirnbaumModel, modelNames=c('Chen', 'Birnbaum'))
genefuModelComparisons <- compareModels(genefuModelComparisons,
    validatedHaiderModel, model2Name='Haider')

allModelComparisons <- compareModels(genefuModelComparisons, 
    validatedPCOSPmodel, model2Name='PCOSP')
# We are only interested in comparing the summaries, so we subset 
#>our model comparison
allModelComparisons <- subset(allModelComparisons, isSummary == TRUE)

allDindexComparisonForestPlot <- forestPlot(allModelComparisons,
    stat='log_D_index', colourBy='model', groupBy='mDataType')

ggsave(allDindexComparisonForestPlot, width=8.5, height=11,
    file='results/figures/classifierlogDindexForestPlot.pdf')

allConcIndexComparisonForestPlot <- forestPlot(allModelComparisons,
    stat='concordance_index', colourBy='model', groupBy='mDataType')

ggsave(allConcIndexComparisonForestPlot, width=8.5, height=11,
    file='results/figures/classifierCindexForestPlot.pdf')

## ---- 6. Comparing PCOSP Models By Patient Subtype
data(cohortSubtypeDFs)

# Add the subtypes to the prediction cohort
subtypedPCOSPValCohorts <- assignSubtypes(PCOSPpredValCohorts, 
    cohortSubtypeDFs)

subtypeValidatedPCOSPmodel <- validateModel(trainedPCOSPmodel, 
    valData=subtypedPCOSPValCohorts)

qsave(subtypeValidatedPCOSPmodel, file='results/validatedSubtype.qs')

subtypeDindexForestPlot <- 
    forestPlot(subtypeValidatedPCOSPmodel, stat='log_D_index', 
        groupBy='cohort', colourBy='subtype')

ggsave(subtypeDindexForestPlot, width=8.5, height=11,
    file='results/figures/subtypeLogDindexForestPlot.pdf')

subtypeCindexForestPlot <- 
    forestPlot(subtypeValidatedPCOSPmodel, stat='concordance_index',
        groupBy='cohort', colourBy='subtype')

ggsave(subtypeCindexForestPlot, width=8.5, height=11,
    file='results/figures/subtypeCindexForestPlot.pdf')
