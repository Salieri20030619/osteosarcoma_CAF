cat('R=',as.character(getRversion()),'\n')
for (p in c('scDblFinder','SingleCellExperiment','zellkonverter')) cat(p,requireNamespace(p,quietly=TRUE),'\n')
