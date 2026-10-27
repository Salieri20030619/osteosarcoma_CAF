# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' install_copykat.R
# RUN_COMMAND: Rscript --vanilla install_copykat.R
# EXPECTED_SUCCESS_MARKER: CNV_ENVIRONMENT_READY
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914'
dir.create(file.path(root,'Rlib'),recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(root,'Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=300,Ncpus=4)
needed <- c('dlm','parallelDist','mixtools','MCMCpack','transport','jsonlite','data.table')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,lib=file.path(root,'Rlib'))
stopifnot(all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))
if(!requireNamespace('copykat',quietly=TRUE)) install.packages(file.path(root,'resources','copykat.tar.gz'),repos=NULL,type='source',lib=file.path(root,'Rlib'))
stopifnot(requireNamespace('copykat',quietly=TRUE),requireNamespace('infercnv',quietly=TRUE))
sink(file.path(root,'logs','cnv_sessionInfo.txt'));print(sessionInfo());print(.libPaths());print(formals(copykat::copykat));sink()
suppressPackageStartupMessages(library(copykat))
full.anno <- get('full.anno',envir=as.environment('package:copykat'))
write.table(full.anno,file.path(root,'resources','copykat_full_anno.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('CNV_ENVIRONMENT_READY\n')
