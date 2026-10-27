O <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase36_external_validation_20260927'
s <- read.delim(file.path(O,'tables/site_comparison_all_prespecified_features.tsv'),check.names=FALSE)
v <- read.delim(file.path(O,'tables/patient_site_feature_values.tsv'),check.names=FALSE)
r <- read.delim(file.path(O,'tables/site_comparison_patient_roster.tsv'),check.names=FALSE)
stopifnot(nrow(s)==130, length(unique(r$patient_id))==15, sum(r$role=='paired')==4, sum(r$role=='independent_primary')==5, sum(r$role=='independent_lung')==6)
out <- s; out$P_independent_R <- NA_real_
for(i in seq_len(nrow(s))){
 z <- v[v$feature==s$feature[i] & v$feature_type==s$feature_type[i],]
 if(s$contrast[i]=='paired_lung_minus_primary'){
  ids <- r$patient_id[r$role=='paired']
  a <- z[z$site_group=='primary_tumor',]; b <- z[z$site_group=='lung_metastasis',]
  d <- b$value[match(ids,b$patient_id)]-a$value[match(ids,a$patient_id)];d<-d[is.finite(d)]
  if(length(d)>=2 && any(d!=0)){
   grid<-as.matrix(expand.grid(rep(list(c(-1,1)),length(d))))
   out$P_independent_R[i]<-mean(abs(as.vector(grid %*% d)/length(d))>=abs(mean(d))-1e-12)
  }
 } else {
  a<-z$value[z$site_group=='primary_tumor' & z$patient_id %in% r$patient_id[r$role=='independent_primary']]
  b<-z$value[z$site_group=='lung_metastasis' & z$patient_id %in% r$patient_id[r$role=='independent_lung']]
  a<-a[is.finite(a)];b<-b[is.finite(b)];all<-c(a,b)
  if(length(a)>=2 && length(b)>=2 && length(unique(all))>1){
   delta<-combn(seq_along(all),length(a),FUN=function(ix) mean(all[-ix])-mean(all[ix]))
   out$P_independent_R[i]<-mean(abs(delta)>=abs(mean(b)-mean(a))-1e-12)
  }
 }
}
stopifnot(identical(is.na(out$P),is.na(out$P_independent_R)))
stopifnot(max(abs(out$P-out$P_independent_R),na.rm=TRUE)<1e-12)
out$FDR_independent_R<-NA_real_
for(fam in unique(out$family)){
 ix<-which(out$family==fam);p<-out$P_independent_R[ix]
 if(any(is.finite(p)))out$FDR_independent_R[ix]<-p.adjust(p,method='BH',n=sum(is.finite(p)))
}
stopifnot(identical(is.na(out$FDR),is.na(out$FDR_independent_R)),max(abs(out$FDR-out$FDR_independent_R),na.rm=TRUE)<1e-12)
write.table(out,file.path(O,'tables/independent_R_site_statistics_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
cat('INDEPENDENT_R_EXACT_P_AND_BH_ACCEPTED',nrow(out),'TEST_ROWS\n')
