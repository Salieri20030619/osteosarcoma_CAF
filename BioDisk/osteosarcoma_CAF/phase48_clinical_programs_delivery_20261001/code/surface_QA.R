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

Q <- file.path(P,'phase48_clinical_programs_delivery_20261001')
qa[[1]] <- checkset(readtsv(file.path(Q,'tables/surface_competing_lineage_statistics.tsv')),readtsv(file.path(Q,'tables/surface_competing_lineage_patient_effects.tsv')),c('layer','class_id','dataset','receiver','gene'),'surface_competing_lineage')
qa[[2]] <- checkset(readtsv(file.path(Q,'tables/external_surface_gene_statistics.tsv')),readtsv(file.path(Q,'tables/external_surface_patient_effects.tsv')),c('layer','class_id','tissue','gene'),'external_surface_genes')
for(name in c('new_surface_pair_RNA_codetection.tsv','external_new_surface_pair_RNA_codetection.tsv')) {
 t<-readtsv(file.path(Q,'tables',name));ok<-is.finite(t$both_detected);stopifnot(all(t$both_detected[ok]>=0),all(t$both_detected[ok]<=t$cells[ok]));nz<-ok&t$cells>0;stopifnot(max(abs(t$fraction_both[nz]-t$both_detected[nz]/t$cells[nz]))<1e-12,all(is.na(t$fraction_both[t$cells==0])))
}
write.table(do.call(rbind,qa),file.path(Q,'tables/independent_R_surface_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(Q,'SURFACE_INDEPENDENT_R_QA_COMPLETE'));print(do.call(rbind,qa))
