# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' install_scevan.R
# RUN_COMMAND: Rscript --vanilla install_scevan.R
# EXPECTED_SUCCESS_MARKER: SCEVAN_ENV_READY
root<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914'
lib<-file.path(root,'SCEVAN_Rlib');dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,file.path(root,'Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600,Ncpus=4)
required<-c('parallelDist','pheatmap','forcats','dplyr','fgsea','cluster','ggplot2','Rtsne','scran','ape','ggtree','tidytree','ggrepel')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
cat('MISSING',paste(missing,collapse=','),'\n')
if(length(missing)) {
 if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',lib=lib)
 BiocManager::install(missing,lib=lib,ask=FALSE,update=FALSE,Ncpus=4)
}
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
if(!requireNamespace('yaGST',quietly=TRUE)) install.packages(file.path(root,'resources','yaGST.tar.gz'),repos=NULL,type='source',lib=lib)
if(!requireNamespace('SCEVAN',quietly=TRUE)) install.packages(file.path(root,'resources','SCEVAN.tar.gz'),repos=NULL,type='source',lib=lib)
suppressPackageStartupMessages(library(SCEVAN))
sink(file.path(root,'logs/scevan_sessionInfo.txt'));print(sessionInfo());print(formals(SCEVAN::pipelineCNA));sink()
cat('SCEVAN_ENV_READY\n')
