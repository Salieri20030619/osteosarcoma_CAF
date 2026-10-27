options(stringsAsFactors=FALSE)
o <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase60_F3_review_20261005'
rt <- function(p) read.delim(p,check.names=FALSE,na.strings=c('','NA','nan'))
verify <- function(kind) {
 s<-rt(file.path(o,'tables',paste0(kind,'_statistics.tsv')));e<-rt(file.path(o,'tables',paste0('patient_',kind,'_effects.tsv')));item<-if(kind=='program')'program' else 'gene'
 key<-function(x)do.call(paste,c(x[c('dataset','matrix','context',item)],sep='|'));sk<-key(s);ek<-key(e);stopifnot(!anyDuplicated(sk),!anyDuplicated(paste(ek,e$donor_record)));ep<-split(e$effect,ek);pe<-c();fe<-c()
 for(i in seq_len(nrow(s))) {
  v<-ep[[sk[i]]];v<-v[is.finite(v)];n<-length(v);stopifnot(n==s$N[i]);if(n)stopifnot(abs(mean(v)-s$mean_delta[i])<1e-9);stopifnot(sum(v>0)==s$positive_patients[i],sum(v<0)==s$negative_patients[i])
  if(n<2 || all(v==0)) {stopifnot(is.na(s$P[i]));next}
  signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));p<-mean(abs(signs%*%v/n)>=abs(mean(v))-1e-12);pe<-c(pe,abs(p-s$P[i]))
 }
 for(f in unique(s$family)) {ii<-which(s$family==f);p<-s$P[ii];p[is.na(p)]<-1;q<-p.adjust(p,'BH');ok<-!is.na(s$P[ii]);fe<-c(fe,abs(q[ok]-s$FDR[ii][ok]));stopifnot(all(is.na(s$FDR[ii][!ok])))}
 stopifnot(max(pe,0)<1e-9,max(fe,0)<1e-9);data.frame(type=kind,rows=nrow(s),valid_P=sum(!is.na(s$P)),families=length(unique(s$family)),max_P_error=max(pe,0),max_FDR_error=max(fe,0))
}
z<-rbind(verify('program'),verify('marker'));write.table(z,file.path(o,'QA/independent_R_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(as.data.frame(z));writeLines('SUCCESS',file.path(o,'QA/R_COMPLETE'))
