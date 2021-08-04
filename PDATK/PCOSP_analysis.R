library(MetaGxPancreas)
library(PDATK)
library(BiocParallel)

# NOTE: This will download several GiB of data to your ExperimentHub Cache
metaGxPanc <- loadPancreasDatasets()
pancreasCohorts <- metaGxPanc$SummarizedExperiments

cohortsWithSurvival <- c("ICGCMICRO", "ICGCSEQ", "PCSI", "TCGA", "KIRBY", "UNC", "CHEN", 
"COLLISON", "ZHANG", "OUH", "WINTER")
names(pancreasCohorts) <- gsub('_SumExp', '', names(pancreasCohorts))
pancreasCohortsWithSurv <- pancreasCohorts[cohortsWithSurvival]

# SurvivalExperiment is just a wrapper around SummarizedExpeirment
pancreasSurvExps <- lapply(pancreasCohortsWithSurv, FUN=SurvivalExperiment, 
    survival_time='days_to_death', event_occurred='vital_status')
# 
pancreasCohortList <- CohortList(pancreasSurvExps)

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

# Set a seed for a reproducible analysis!
randomSeed <- 1987
set.seed(randomSeed)

# Construct a PCOSP model object
PCOSPmodel <- PCOSP(ICGCtrainCohorts, minDaySurvived=365,
    randomSeed=randomSeed)

trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=15, minAccuracy=0.6)

metadata(trainedPCOSPmodel)$modelParams