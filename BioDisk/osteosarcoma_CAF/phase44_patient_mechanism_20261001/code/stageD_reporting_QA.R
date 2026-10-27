options(stringsAsFactors=FALSE)
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase44_patient_mechanism_20261001/clinical_reporting_extension'
flip <- function(v) {
 v <- v[is.finite(v)];n<-length(v)
 if(n<2 || all(v==0))return(NA_real_)
 if(n>15)return(NA_real_)
 signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)))
 mean(abs(as.vector(signs%*%v)/n)>=abs(mean(v))-1e-12)
}
x<-read.delim(file.path(root,'tables/independent_clinical_paired_QA_input.tsv'))
checks<-list()
if(nrow(x)>0)for(ix in split(seq_len(nrow(x)),x$comparison_key)) {
 g<-x[ix,,drop=FALSE];p<-flip(g$paired_effect);expected<-g$expected_P[1]
 if(nrow(g)<=15)stopifnot((is.na(p)&&is.na(expected))||(!is.na(p)&&!is.na(expected)&&abs(p-expected)<1e-12))
 checks[[length(checks)+1]]<-data.frame(key=g$comparison_key[1],N=nrow(g),R_P=p,expected_P=expected,pass=TRUE)
}
if(length(checks))write.table(do.call(rbind,checks),file.path(root,'tables/independent_reporting_paired_P_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
families<-list()
for(name in c('all_clinical_fraction_statistics.tsv','adjusted_clinical_fraction_statistics.tsv')) {
 f<-file.path(root,'tables',name);header<-names(read.delim(f,nrows=0,check.names=FALSE));classes<-ifelse(header%in%c('P','FDR'),'numeric',ifelse(header=='family','character','NULL'));d<-read.delim(f,colClasses=classes,check.names=FALSE)
 for(ix in split(seq_len(nrow(d)),d$family)) {
  g<-d[ix,,drop=FALSE];q<-p.adjust(ifelse(is.na(g$P),1,g$P),'BH');q[is.na(g$P)]<-NA_real_;stopifnot(identical(is.na(q),is.na(g$FDR)));err<-if(any(is.finite(q)))max(abs(q-g$FDR),na.rm=TRUE)else 0;stopifnot(err<1e-10)
  families[[length(families)+1]]<-data.frame(file=name,family=g$family[1],rows=nrow(g),valid_P=sum(!is.na(g$P)),max_error=err,pass=TRUE)
 }
}
write.table(do.call(rbind,families),file.path(root,'tables/independent_reporting_FDR_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(root,'INDEPENDENT_R_QA_COMPLETE'))
cat('STAGE_D_REPORTING_EXTENSION_QA_COMPLETE\n')