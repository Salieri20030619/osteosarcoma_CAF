suppressPackageStartupMessages(library(infercnv))
ns<-asNamespace('infercnv')
for(n in ls(ns,all.names=TRUE)[grepl('meanvar|mean_var',ls(ns,all.names=TRUE))]) {cat('FUNCTION',n,'\n');cat(paste(deparse(get(n,envir=ns)),collapse='\n'),'\n')}
