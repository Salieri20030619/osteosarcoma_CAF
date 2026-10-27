
.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
for(pkg in c('copykat','infercnv')) {
 ns<-asNamespace(pkg)
 for(nm in ls(ns,all.names=TRUE)) {
  f<-get(nm,envir=ns)
  if(!is.function(f))next
  lines<-deparse(f)
  hits<-grep('mclapply|makeCluster|parLapply|num_threads|n.cores|registerDoParallel',lines)
  if(length(hits)) {
   cat('\n',pkg,'::',nm,'\n',sep='')
   cat(paste(lines[hits],collapse='\n'),'\n')
  }
 }
}
