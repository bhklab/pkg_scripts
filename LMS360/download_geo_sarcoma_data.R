library(GEOquery)
library(Biobase)
library(SummarizedExperiment)

# install appropriate brain array annotations for array platform
brain_array_urls <- function(array, species="hs", annotation="ensg", version="25.0.0") {
    ## FIXME:: make robust to missing arguments
    paste0("http://mbni.org/customcdf/", version, "/", annotation, ".download",
        "/", c("", "", "pd."), array, c("", "", "."), species, c("", "", "."),
        annotation, c("cdf", "probe", ""), "_", version, ".tar.gz")
}
array <- "hgu133a"
brain_array <- brain_array_urls(array)
for (pkg in brain_array) {
    install.packages(pkg, type="src", repos=NULL)
}

datasets <- c("GSE21122", "GSE21050")