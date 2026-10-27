options(stringsAsFactors=FALSE)
P <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
B <- file.path(P,'phase56_patient_programs_20261004')
O <- file.path(P,'phase57_downstream_update_20261004')
Q <- file.path(O,'independent_QA')
rd <- function(f) read.delim(if(grepl('gz$',f)) gzfile(f) else f,check.names=FALSE)
cache <- new.env()
signs_for <- function(n) {
 k<-as.character(n)
 if(exists(k,envir=cache,inherits=FALSE))return(get(k,envir=cache))
 if(n<=15) { ids<-0:(2^n-1); s<-sapply(0:(n-1),function(j) ifelse(bitwAnd(as.integer(ids),bitwShiftL(1L,j))>0,1,-1)) }
 else s<-as.matrix(read.table(gzfile(file.path(Q,paste0('MC_signs_N',n,'.tsv.gz'))),header=FALSE))
 assign(k,s,envir=cache);s
}
calcP <- function(d) {
 out<-rep(NA_real_,ncol(d));n<-nrow(d)
 if(n<2)return(out)
 keep<-colSums(is.finite(d))==n & colSums(abs(d),na.rm=TRUE)>0
 if(!any(keep))return(out)
 z<-d[,keep,drop=FALSE];ss<-signs_for(n);ob<-abs(colMeans(z));null<-abs(ss %*% z/n)
 hits<-colSums(sweep(null,2,ob-1e-12,'>='));out[keep]<-if(n<=15)hits/nrow(ss) else (hits+1)/(nrow(ss)+1);out
}
check <- function(stats,eff,itemcol,label) {
 result<-list();i<-0
 for(ctx in unique(stats$context_id)) {
  g<-stats[stats$context_id==ctx,,drop=FALSE];e<-eff[eff$context_id==ctx,,drop=FALSE];donors<-sort(unique(e$donor_record));items<-g[[itemcol]]
  stopifnot(!anyDuplicated(items),!anyDuplicated(paste(e$donor_record,e[[itemcol]],sep='|')))
  d<-matrix(NA_real_,length(donors),length(items));d[cbind(match(e$donor_record,donors),match(e[[itemcol]],items))]<-e$effect
  pp<-calcP(d);stopifnot(identical(is.na(pp),is.na(g$P)))
  er<-if(any(is.finite(pp)))max(abs(pp-g$P),na.rm=TRUE) else 0
  stopifnot(er<1e-10)
  i<-i+1;result[[i]]<-data.frame(type=label,context_id=ctx,N=length(donors),rows=nrow(g),max_P_error=er)
 }
 do.call(rbind,result)
}
s<-rd(file.path(B,'tables/all_patient_program_statistics.tsv.gz'));e<-rd(file.path(B,'tables/patient_program_effects.tsv.gz'))
a<-check(s,e,'program','all_program_P')
g<-rd(file.path(Q,'gene_P_check_statistics.tsv'));ge<-rd(file.path(Q,'gene_P_check_effects.tsv'));b<-check(g,ge,'gene','selected_gene_P')
write.table(rbind(a,b),file.path(Q,'independent_R_P_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
fdrcheck<-function(x,id) {
 expected<-p.adjust(ifelse(is.na(x$P),1,x$P),method='BH');expected[is.na(x$P)]<-NA_real_
 stopifnot(identical(is.na(expected),is.na(x$FDR)))
 err<-if(any(is.finite(expected)))max(abs(expected-x$FDR),na.rm=TRUE) else 0
 stopifnot(err<1e-10);data.frame(family=id,rows=nrow(x),max_FDR_error=err)
}
fdr<-lapply(split(s,s$family),function(x)fdrcheck(x,x$family[1]))
for(f in c(list.files(file.path(O,'tables/genes'),pattern='tsv.gz$',full.names=TRUE),list.files(file.path(O,'tables/pathways'),pattern='tsv.gz$',full.names=TRUE))) {
 x<-rd(f);fdr[[length(fdr)+1]]<-fdrcheck(x,x$family[1])
}
write.table(do.call(rbind,fdr),file.path(Q,'independent_R_FDR_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('SUCCESS','All patient-program P values, selected gene P values, all gene/program/ORA FDR families checked in R.'),file.path(Q,'INDEPENDENT_R_QA_COMPLETE'))
cat('INDEPENDENT_R_QA_COMPLETE',nrow(a),nrow(b),length(fdr),'\n')