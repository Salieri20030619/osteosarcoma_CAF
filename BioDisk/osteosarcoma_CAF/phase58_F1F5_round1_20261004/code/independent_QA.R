options(stringsAsFactors=FALSE)
p<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase58_F1F5_round1_20261004';q<-file.path(p,'independent_QA')
rd<-function(f)read.delim(if(grepl('gz$',f))gzfile(f) else f,check.names=FALSE)
ex<-rd(file.path(q,'P_expected.tsv'));eff<-rd(file.path(q,'P_patient_effects.tsv.gz'));cache<-new.env()
sgn<-function(n){k<-as.character(n);if(exists(k,envir=cache,inherits=FALSE))return(get(k,envir=cache));s<-if(n<=15)sapply(0:(n-1),function(j)ifelse(bitwAnd(as.integer(0:(2^n-1)),bitwShiftL(1L,j))>0,1,-1)) else as.matrix(read.table(gzfile(file.path(q,paste0('MC_signs_N',n,'.tsv.gz'))),header=FALSE));assign(k,s,envir=cache);s}
calc<-function(v){v<-v[is.finite(v)];n<-length(v);if(n<2 || !any(v!=0))return(NA_real_);ss<-sgn(n);h<-sum(abs(ss%*%v/n)>=abs(mean(v))-1e-12);if(n<=15)h/nrow(ss) else (h+1)/(nrow(ss)+1)}
res<-list();idx<-0
byctx<-split(eff,eff$context_id)
for(ctx in unique(ex$context_id)) {
 z<-ex[ex$context_id==ctx,,drop=FALSE];e<-byctx[[ctx]]
 byitem<-if(is.null(e))list() else split(e$effect,e$item)
 pp<-vapply(z$item,function(item)calc(byitem[[item]]),numeric(1));nn<-vapply(z$item,function(item)sum(is.finite(byitem[[item]])),integer(1))
 stopifnot(all(nn==z$N),all(is.na(pp)==is.na(z$P)));err<-if(any(is.finite(pp)))max(abs(pp-z$P),na.rm=TRUE) else 0;stopifnot(err<1e-10)
 idx<-idx+1;res[[idx]]<-data.frame(context_id=ctx,rows=nrow(z),max_P_error=err)
}
write.table(do.call(rbind,res),file.path(q,'independent_R_P_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
fd<-list();add<-function(x,label,pcol='P',qcol='FDR'){
 pv<-x[[pcol]];expected<-p.adjust(ifelse(is.na(pv),1,pv),method='BH');expected[is.na(pv)]<-NA_real_;actual<-x[[qcol]];stopifnot(identical(is.na(expected),is.na(actual)));err<-if(any(is.finite(expected)))max(abs(expected-actual),na.rm=TRUE) else 0;stopifnot(err<1e-10);data.frame(family=label,rows=nrow(x),max_FDR_error=err)
}
files<-c(list.files(file.path(p,'tables/genes'),pattern='tsv.gz$',full.names=TRUE),file.path(p,'tables/all_program_statistics.tsv.gz'),file.path(p,'tables/corrected_program_statistics.tsv.gz'),file.path(p,'tables/surface_competition_statistics.tsv'),list.files(file.path(p,'tables/ora'),pattern='tsv.gz$',full.names=TRUE))
for(f in files){x<-rd(f);for(g in split(x,x$family))fd[[length(fd)+1]]<-add(g,g$family[1]);rm(x);gc(verbose=FALSE)}
for(f in list.files(file.path(p,'tables/gsea'),pattern='tsv.gz$',full.names=TRUE)){x<-rd(f);fd[[length(fd)+1]]<-add(x,basename(f),'GSEA_model_P','model_BH_FDR')}
write.table(do.call(rbind,fd),file.path(q,'independent_R_FDR_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines('SUCCESS',file.path(q,'INDEPENDENT_R_COMPLETE'));cat('INDEPENDENT_R_COMPLETE',nrow(ex),length(fd),'\n')