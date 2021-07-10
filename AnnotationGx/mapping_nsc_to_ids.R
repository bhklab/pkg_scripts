library(AnnotationGx)
library(data.table)

# -- mapping NSC numbers to CIDs and drug names
NCI60 <- fread('local_data/DTP_NCI60_RAW.csv')
NCI60[`PubChem SID` == -1, `PubChem SID` := NA]
ids <- unique(na.omit(NCI60[[1]]))
NSCtoCID <- getPubChemFromNSC(ids[1:10], proxy=FALSE)

# -- retry batch queries until no more matches are returned
failed <- getFailed(NSCtoCID)
failedIDs <- getFailedIDs(NSCtoCID)
if (length(failedIDs > 1)) {
    retryQueries <- getPubChemFromNSC(failedIDs)
    NSCtoCID <- rbind(NSCtoCID, retryQueries)
    while (nrow(retryQueries) > 0) {
        failedIDs <- getFailedIDs(retryQueries)
        retryQueries <- getPubChemFromNSC(failedIDs)
        NSCtoCID <- rbind(NSCtoCID, retryQueries)
    }
}
# -- retry non-batch queries until no matches are returned
failedIDs <- getFailedIDs(retryQueries)
failedMsg <- getFailureMessages(retryQueries)
if (all.equal(failedMsg[, unique(Code)], 'PUGREST.BadRequest'))
if (length(failedIDs > 1)) {
    retryQueries2 <- getPubChemFromNSC(failedIDs, batch=FALSE)
    failedIDs <- getFailedIDs(retryQueries)
    while (nrow(retryQueries2) > 0) {
        retryQueries2 <- getPubChemFromNSC(failedIDs, batch=FALSE)
        failedIDs <- getFailedIDs(retryQueries2)
        NSCtoCID <- rbind(NSCtoCID, retryQueries2)
    }
}
failedIDs <- getFailedIDs(retryQueries2)

# -- try to 
colnames(NCI60) <- gsub('_#', '', gsub(' ', '_', colnames(NCI60)))
oldSIDs <- 
cids <- na.omit(NSCtoCID$CID)
compoundProperties <- getPubChemCompound(cids)
unmapped_sids <- NSCtoCID[!(CID %chin% compoundProperties$CID), ]$SID
compoundPropertiesSID <- getPubChemCompound(unmapped_sids, from='sid')
NSCtoName <- merge.data.table(NSCtoCID, compoundProperties, by='CID', all.x=TRUE)
NSCtoNameSID <- merge.data.table(NSCtoCID, compoundPropertiesSID, by='SID')
NSCtoName[NSCtoNameSID, Title := i.Title, on=c('SID', 'NSC_id')]
failed_to_map <- setdiff(ids, NSCtoName$NSC_id)
retryAgain <- getPubChemFromNSC(failed_to_map, to='sids', batch=FALSE)
SIDtoName <- getPubChemCompound(sids, from='sid')