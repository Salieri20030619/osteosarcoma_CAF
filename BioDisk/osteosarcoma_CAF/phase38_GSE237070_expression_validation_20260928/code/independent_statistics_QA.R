O <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase38_GSE237070_expression_validation_20260928"
t <- read.delim(file.path(O,"tables/patient_level_marker_program_statistics.tsv"), check.names=FALSE)
e <- read.delim(file.path(O,"tables/individual_patient_effects.tsv"),check.names=FALSE)
k <- c("grouping","tissue","left","right","kind","endpoint")
key <- function(x) do.call(paste,c(x[k],sep="|"))
te <- key(t);ee<-key(e);p<-rep(NA_real_,nrow(t))
for(i in seq_len(nrow(t))) {
 v<-e$effect[ee==te[i]]
 stopifnot(length(v)==t$N_patients[i])
 if(length(v)>=2 && !is.na(t$P[i])) {
  signs<-as.matrix(expand.grid(rep(list(c(-1,1)),length(v))))
  stat<-abs(mean(v));null<-abs(as.vector(signs %*% v)/length(v))
  p[i]<-mean(null>=stat-1e-12)
 }
}
stopifnot(identical(is.na(p),is.na(t$P)),max(abs(p-t$P),na.rm=TRUE)<1e-12)
q<-rep(NA_real_,nrow(t))
for(f in unique(t$family)) {
 ix<-which(t$family==f);x<-p[ix];x[is.na(x)]<-1;qq<-p.adjust(x,"BH");qq[is.na(p[ix])]<-NA_real_;q[ix]<-qq
}
stopifnot(identical(is.na(q),is.na(t$FDR)),max(abs(q-t$FDR),na.rm=TRUE)<1e-12)
write.table(data.frame(rows=nrow(t),valid_tests=sum(!is.na(p)),P_max_abs_difference=max(abs(p-t$P),na.rm=TRUE),FDR_max_abs_difference=max(abs(q-t$FDR),na.rm=TRUE),native_R_QA="PASS"),file.path(O,"tables/independent_R_statistics_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
