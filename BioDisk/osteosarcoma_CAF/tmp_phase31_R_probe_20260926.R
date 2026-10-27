 .libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
for(pkg in c('decontX','celda','SingleCellExperiment','scDblFinder','Seurat','Matrix','infercnv','jsonlite')) cat(pkg,requireNamespace(pkg,quietly=TRUE),if(requireNamespace(pkg,quietly=TRUE))as.character(packageVersion(pkg)) else '', '\n')
if(requireNamespace('decontX',quietly=TRUE))print(formals(decontX::decontX))
