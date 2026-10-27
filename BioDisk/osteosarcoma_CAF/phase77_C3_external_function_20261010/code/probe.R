cat(R.version.string,"\n"); for(p in c("qs2","SeuratObject","Matrix")) cat(p,requireNamespace(p,quietly=TRUE),"\n")
