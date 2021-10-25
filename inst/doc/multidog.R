## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=4.5,
  fig.height=3.5
)

## ----setup--------------------------------------------------------------------
library(future)
library(updog)
data("uitdewilligen")

## -----------------------------------------------------------------------------
refmat  <- t(uitdewilligen$refmat)
sizemat <- t(uitdewilligen$sizemat)
ploidy  <- uitdewilligen$ploidy

## -----------------------------------------------------------------------------
setdiff(colnames(sizemat), colnames(refmat))
setdiff(rownames(sizemat), rownames(refmat))

## -----------------------------------------------------------------------------
future::availableCores()

## -----------------------------------------------------------------------------
mout <- multidog(refmat = refmat, 
                 sizemat = sizemat, 
                 ploidy = ploidy, 
                 model = "norm",
                 nc = 2)

## ---- eval = FALSE------------------------------------------------------------
#  future::plan(future::multisession, workers = nc)

## ---- eval = FALSE------------------------------------------------------------
#  future::plan(future::multicore, workers = 2)
#  mout <- multidog(refmat = refmat,
#                   sizemat = sizemat,
#                   ploidy = ploidy,
#                   model = "norm",
#                   nc = NA)
#  
#  ## Shut down parallel workers
#  future::plan(future::sequential)

## -----------------------------------------------------------------------------
plot(mout, indices = c(1, 5, 100))

## -----------------------------------------------------------------------------
str(mout$snpdf)

## -----------------------------------------------------------------------------
str(mout$inddf)

## -----------------------------------------------------------------------------
genomat <- format_multidog(mout, varname = "geno")
head(genomat)

## -----------------------------------------------------------------------------
dim(mout$snpdf)
dim(mout$inddf)
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > exp(-1) & bias < exp(1))
dim(mout_cleaned$snpdf)
dim(mout_cleaned$inddf)

## ---- eval = FALSE------------------------------------------------------------
#  # install.packages("BiocManager")
#  # BiocManager::install(c("VariantAnnotation", "GenomicRanges", "S4Vectors", "IRanges"))
#  export_vcf(obj = mout, filename = "./multidog_fit.vcf")

