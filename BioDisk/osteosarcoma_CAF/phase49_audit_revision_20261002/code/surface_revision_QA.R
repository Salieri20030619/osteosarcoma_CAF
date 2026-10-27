P <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase49_audit_revision_20261002'
S <- file.path(P,'surface','tables')
readtab <- function(f) read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
check <- function(st,pa,keys,label) {
 vectors <- lapply(seq_len(nrow(st)), function(i) { ok <- rep(TRUE,nrow(pa)); for(k in keys) ok <- ok & pa[[k]]==st[[k]][i]; v<-pa$effect[ok];v[is.finite(v)] })
 pp <- rep(NA_real_,nrow(st));meanv<-vapply(vectors,function(v)if(length(v))mean(v)else NA_real_,numeric(1))
 for(n in unique(lengths(vectors))) {
  ids<-which(lengths(vectors)==n & vapply(vectors,function(v)any(v!=0),logical(1)));if(n<2 || !length(ids))next;stopifnot(n<=20)
  mat<-do.call(cbind,vectors[ids]);obs<-abs(colMeans(mat));hits<-numeric(length(ids))
  for(start in seq(0,2^n-1,by=4096)) { ints<-as.integer(start:min(start+4095,2^n-1));signs<-sapply(seq_len(n),function(j)ifelse(bitwAnd(ints,bitwShiftL(1L,j-1L))!=0,1,-1));null<-abs(signs%*%mat/n);hits<-hits+colSums(sweep(null,2,obs-1e-12,FUN='>=')) }
  pp[ids]<-hits/2^n
 }
 stopifnot(identical(is.na(pp),is.na(st$P)),all(lengths(vectors)==st$N));q<-rep(NA_real_,nrow(st))
 for(f in unique(st$family)) { ix<-which(st$family==f);v<-pp[ix];ok<-is.finite(v);v[!ok]<-1;qq<-p.adjust(v,'BH');qq[!ok]<-NA;q[ix]<-qq }
 errs<-c(max(c(0,abs(pp-st$P)),na.rm=TRUE),max(c(0,abs(q-st$FDR)),na.rm=TRUE),max(c(0,abs(meanv-st$mean_effect)),na.rm=TRUE));stopifnot(all(errs<1e-10));data.frame(module=label,rows=nrow(st),P_error=errs[1],FDR_error=errs[2],mean_error=errs[3],passed=TRUE)
}
a<-check(readtab(file.path(S,'surface_competing_lineage_statistics.tsv')),readtab(file.path(S,'surface_competing_lineage_patient_effects.tsv')),c('layer','class_id','dataset','receiver','gene'),'surface_competing')
b<-check(readtab(file.path(S,'external_surface_gene_statistics.tsv')),readtab(file.path(S,'external_surface_patient_effects.tsv')),c('layer','class_id','tissue','gene'),'surface_external')
co<-readtab(file.path(S,'new_surface_pair_RNA_codetection.tsv'));stopifnot(all(co$both_detected[is.finite(co$both_detected)]>=0),all(co$both_detected[is.finite(co$both_detected)]<=co$cells[is.finite(co$both_detected)]),length(unique(co$pair))==45)
write.table(rbind(a,b),file.path(P,'surface','delivery','INDEPENDENT_R_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE);writeLines('SUCCESS',file.path(P,'surface','SURFACE_INDEPENDENT_R_QA_COMPLETE'));print(as.data.frame(rbind(a,b)))
