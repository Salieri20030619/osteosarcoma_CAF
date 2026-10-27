lib <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib"
.libPaths(c(lib,.libPaths()))
base <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/R_sources"
for (f in c("RcppParallel_6.2.1.tar.gz","stringfish_0.19.2.tar.gz","qs2_0.3.1.tar.gz")) install.packages(file.path(base,f),repos=NULL,type="source",lib=lib,Ncpus=2)
stopifnot(requireNamespace("qs2",quietly=TRUE))
cat("QS2_READY",as.character(packageVersion("qs2")),"\n")
