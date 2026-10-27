.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
for(p in c('Matrix','data.table','jsonlite','copykat','infercnv')){stopifnot(requireNamespace(p,quietly=TRUE));cat(p,as.character(packageVersion(p)),'\n')}
invisible(parse(file='/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase26_client_reanalysis_20260922/code/native_cnv.R'))
cat('R_PREFLIGHT_OK\n')
