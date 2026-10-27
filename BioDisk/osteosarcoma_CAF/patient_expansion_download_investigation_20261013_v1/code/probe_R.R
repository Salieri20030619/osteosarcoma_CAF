cat(R.version.string,"\n")
for (p in c("qs2","Seurat","SeuratObject","Matrix","jsonlite")) cat(p,requireNamespace(p,quietly=TRUE),"\n")
