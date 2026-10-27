P <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
O <- file.path(P,'phase46_spatial_functional_20261001')
E <- file.path(P,'phase47_external_candidate_integration_20261001')
readtsv <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
qa <- list()
checkset <- function(stats,pat,keys,label) {
 vals <- lapply(seq_len(nrow(stats)),function(i) {
  keep <- rep(TRUE,nrow(pat))
  for(k in keys) keep <- keep & pat[[k]]==stats[[k]][i]
  d <- pat$effect[keep];d[is.finite(d)]
 })
 pp <- rep(NA_real_,nrow(stats))
 for(n in unique(lengths(vals))) {
  ids <- which(lengths(vals)==n & vapply(vals,function(v) any(v!=0),logical(1)))
  if(n<2 || !length(ids)) next
  stopifnot(n<=20)
  mat <- do.call(cbind,vals[ids]);obs <- abs(colMeans(mat));hits <- rep(0,length(ids))
  for(start in seq(0,2^n-1,by=16384)) {
   ints <- as.integer(start:min(start+16383,2^n-1))
   signs <- sapply(seq_len(n),function(j) ifelse(bitwAnd(ints,bitwShiftL(1L,j-1L))!=0,1,-1))
   null <- abs(signs %*% mat/n)
   hits <- hits+colSums(sweep(null,2,obs-1e-12,FUN='>=') )
  }
  pp[ids] <- hits/2^n
 }
 stopifnot(identical(is.na(pp),is.na(stats$P)))
 err <- max(c(0,abs(pp-stats$P)),na.rm=TRUE);stopifnot(err<1e-12)
 stopifnot(all(lengths(vals)==stats$N))
 means <- vapply(vals,function(v) if(length(v)) mean(v) else NA_real_,numeric(1));stopifnot(max(c(0,abs(means-stats$mean_effect)),na.rm=TRUE)<1e-10)
 q <- rep(NA_real_,nrow(stats))
 for(f in unique(stats$family)) {
  ix <- which(stats$family==f);v <- stats$P[ix];ok <- is.finite(v);v[!ok] <- 1;qq <- p.adjust(v,'BH');qq[!ok] <- NA_real_;q[ix] <- qq
 }
 stopifnot(identical(is.na(q),is.na(stats$FDR)));qe <- max(c(0,abs(q-stats$FDR)),na.rm=TRUE);stopifnot(qe<1e-12)
 data.frame(analysis=label,hypotheses=nrow(stats),P_max_error=err,FDR_max_error=qe,passed=TRUE)
}
qa[[1]] <- checkset(readtsv(file.path(O,'tables/spatial_patient_statistics.tsv')),readtsv(file.path(O,'tables/patient_spatial_effects.tsv')),c('program','endpoint'),'spatial_main')
qa[[2]] <- checkset(readtsv(file.path(O,'tables/spatial_site_sensitivity_statistics.tsv')),readtsv(file.path(O,'tables/patient_site_spatial_effects.tsv')),c('site_group','program','endpoint'),'spatial_site')
qa[[3]] <- checkset(readtsv(file.path(O,'tables/spatial_model_sensitivity_statistics.tsv')),readtsv(file.path(O,'tables/spatial_sensitivity_patient_effects.tsv')),c('program','model'),'spatial_model_sensitivity')
qa[[4]] <- checkset(readtsv(file.path(E,'tables/external_patient_statistics.tsv')),readtsv(file.path(E,'tables/external_patient_effects.tsv')),c('tissue','layer','class_id','kind','endpoint'),'external_migration')
qa[[5]] <- checkset(readtsv(file.path(E,'tables/source_projection_statistics.tsv')),readtsv(file.path(E,'tables/source_projection_patient_effects.tsv')),c('layer','dataset','class_id','kind','endpoint'),'source_projection')
f <- gzfile(file.path(O,'tables/all_spot_RNA_and_components.tsv.gz'),open='rt');z <- read.delim(f,check.names=FALSE,stringsAsFactors=FALSE);close(f)
raw <- readtsv(file.path(O,'tables/independent_QA_raw_gene_counts.tsv'));ed <- readtsv(file.path(O,'tables/independent_QA_adjacency.tsv'))
partial <- function(x,y,c) {
 ok <- is.finite(x)&is.finite(y)&apply(is.finite(c),1,all);x<-rank(x[ok]);y<-rank(y[ok]);c<-apply(c[ok,,drop=FALSE],2,rank);m<-cbind(1,c);cor(lm.fit(m,x)$residuals,lm.fit(m,y)$residuals)
}
sample <- readtsv(file.path(O,'tables/sample_spatial_associations.tsv'))
for(sm in unique(raw$sample)) {
 t<-raw[raw$sample==sm,];v<-z[match(t$spot_id,z$spot_id),];stopifnot(identical(t$spot_id,v$spot_id))
 score<-rowMeans(log1p(as.matrix(t[,c('COL1A1','LOX','MMP14')])/t$library_UMI*1e4));cd8<-rowMeans(log1p(as.matrix(t[,c('CD8A','CD8B')])/t$library_UMI*1e4));stopifnot(max(abs(score-v$F3_SELECTED_REMODELING))<1e-10,max(abs(cd8-v$CD8_RNA))<1e-10)
 e<-ed[ed$sample==sm,];neighbor<-rep(NA_real_,nrow(v));means<-tapply(cd8[e$target],e$source,mean);neighbor[as.integer(names(means))]<-means;stopifnot(max(abs(neighbor-v$neighbor_CD8_RNA),na.rm=TRUE)<1e-10)
 controls<-cbind(v$cancer,v$fibroblast,log1p(v$library_UMI),v$x,v$y);p<-partial(score,cd8,controls);expected<-sample$effect[sample$sample==sm&sample$program=='F3_SELECTED_REMODELING'&sample$endpoint=='CD8_adjusted'];stopifnot(abs(p-expected)<1e-10)
}
qa[[6]] <- data.frame(analysis='raw_RNA_neighbor_partial_reconstruction',hypotheses=nrow(raw),P_max_error=NA_real_,FDR_max_error=NA_real_,passed=TRUE)
write.table(do.call(rbind,qa),file.path(O,'tables/independent_R_QA_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(O,'INDEPENDENT_R_QA_COMPLETE'));writeLines('SUCCESS',file.path(E,'INDEPENDENT_R_QA_COMPLETE'));print(do.call(rbind,qa))
