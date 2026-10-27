.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages(library(copykat))
cat('copykat',as.character(packageVersion('copykat')),'R',as.character(getRversion()),'R.home',R.home(),'\n')
invisible(parse(file='/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase92_copykat_validation_20261015/code/copykat_full.R'))
