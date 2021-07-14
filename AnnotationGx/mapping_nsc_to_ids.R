library(AnnotationGx)
library(data.table)

# ---- mapping NSC numbers to CIDs
NCI60 <- fread('local_data/DTP_NCI60_RAW.csv')
NCI60[`PubChem SID` == -1, `PubChem SID` := NA]
ids <- unique(na.omit(NCI60[[1]]))
NSCtoCID <- getPubChemFromNSC(ids, proxy=TRUE)

# -- retry batch queries until no more matches are returned
failed <- getFailed(NSCtoCID)
failedIDs <- getFailedIDs(NSCtoCID)
if (length(failedIDs) > 1) { 
    retryQueries <- getPubChemFromNSC(failedIDs, proxy=TRUE)
    NSCtoCID <- rbind(NSCtoCID, retryQueries)
    while (nrow(retryQueries) > 0) {
        failedIDs <- getFailedIDs(retryQueries)
        retryQueries <- getPubChemFromNSC(failedIDs, proxy=TRUE)
        NSCtoCID <- rbind(NSCtoCID, retryQueries)
    }
    failedIDs <- getFailedIDs(retryQueries)
    failedMsg <- getFailureMessages(retryQueries)
}

# -- retry non-batch queries until no matches are returned
if (length(failedIDs) > 1) {
    retryQueries2 <- getPubChemFromNSC(failedIDs, batch=FALSE, proxy=TRUE)
    failedIDs <- getFailedIDs(retryQueries)
    while (nrow(retryQueries2) > 0) {
        retryQueries2 <- getPubChemFromNSC(failedIDs, batch=FALSE, proxy=TRUE)
        failedIDs <- getFailedIDs(retryQueries2)
        NSCtoCID <- rbind(NSCtoCID, retryQueries2)
    }
    failedIDs <- getFailedIDs(retryQueries2)
}

# ---- map from CID to Title, IUPACName, CanonicalSMILES, IsomericSMILES and InChIKey
cids <- unique(na.omit(NSCtoCID$CID))
CIDtoProperties <- getPubChemCompound(cids, proxy=TRUE,
    properties=c('Title', 'IUPACName', 'CanonicalSMILES', 'IsomericSMILES',
        'InChIKey'))

# -- retry CID to properties until no more results are returned
failedIDs <- getFailedIDs(CIDtoProperties)
.failedCIDs <- failedIDs # keep old version incase we need to restart
if (length(failedIDs) > 0) {
    retryProperties <- getPubChemCompound(failedIDs, proxy=TRUE,
        properties=c('Title', 'IUPACName', 'CanonicalSMILES', 'IsomericSMILES',
             'InChIKey'))
    failedIDs <- getFailedIDs(retryProperties)
    while (nrow(retryProperties) > 0) {
        retryProperties <- getPubChemCompound(failedIDs, proxy=TRUE,
            properties=c('Title', 'IUPACName', 'CanonicalSMILES', 
            'IsomericSMILES', 'InChIKey'))
        CIDtoProperties <- rbind(CIDtoProperties, retryProperties)
        failedIDs <- getFailedIDs(retryProperties)
    }
}

# -- map from SID to properties for 