.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths()))
for (pkg in c('decontX','celda','SingleCellExperiment','Matrix','BiocManager')) cat(pkg,requireNamespace(pkg,quietly=TRUE),'\n')
