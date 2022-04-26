library(PharmacoGx)
library(MultiAssayExperiment)

data(GDSCsmall)
pSet <- GDSCsmall

MAE <- pSet |> molecularProfilesSlot() |> MultiAssayExperiment()
LT <- pSet |>
    CoreGx:::.sensitivityToLongTable() |> # this doesn't work correctly
    as("TreatmentResponseExperiment")

pSet2 <- PharmacoSet2(
    name=name(pSet),
    treatment=drugInfo(pSet),
    sample=sampleInfo(pSet),
    molecularProfiles=MAE,
    treatmentResponse=LT,
    curation=curation(pSet)
)