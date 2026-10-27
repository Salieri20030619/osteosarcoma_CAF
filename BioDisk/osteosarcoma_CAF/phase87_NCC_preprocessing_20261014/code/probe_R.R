for (p in c("scDblFinder","SingleCellExperiment","Matrix","BiocParallel")) { stopifnot(requireNamespace(p,quietly=TRUE)); cat(p,as.character(packageVersion(p)),"\n") }
