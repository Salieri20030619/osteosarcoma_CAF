suppressPackageStartupMessages(library(infercnv))
ns<-asNamespace('infercnv')
nm<-ls(ns,all.names=TRUE)
for(n in nm[grepl('hspike',nm)]) {cat('FUNCTION',n,'\n');cat(paste(deparse(get(n,envir=ns)),collapse='\n'),'\n')}
