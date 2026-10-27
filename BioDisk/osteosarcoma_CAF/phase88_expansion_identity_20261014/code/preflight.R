args <- commandArgs(TRUE)
stopifnot(length(args)==2)
mode <- args[[1]]; out <- args[[2]]
if (mode=='cnv') {
 .libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
 packages <- c('infercnv','copykat','Matrix','data.table','jsonlite')
} else packages <- c('decontX','SingleCellExperiment','scDblFinder','BiocParallel','Matrix','jsonlite')
for (p in packages) {
 suppressPackageStartupMessages(library(p,character.only=TRUE))
 cat(p,as.character(packageVersion(p)),'\n')
}
for (f in c('run_cnv.R','run_decontx.R')) invisible(parse(file=file.path(out,'code',f)))
jsonlite::write_json(list(status='PASS',mode=mode,R=as.character(getRversion()),versions=setNames(lapply(packages,function(p)as.character(packageVersion(p))),packages)),file.path(out,'resources',paste0('PREFLIGHT_',mode,'.json')),pretty=TRUE,auto_unbox=TRUE)
cat('PREFLIGHT_PASS',mode,'\n')
