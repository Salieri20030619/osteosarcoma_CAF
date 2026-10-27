lib <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib"
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=600)
if (!requireNamespace("qs2",quietly=TRUE)) install.packages("qs2",lib=lib,Ncpus=2)
stopifnot(requireNamespace("qs2",quietly=TRUE))
cat("QS2_READY",as.character(packageVersion("qs2")),"\n")
