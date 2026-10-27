options(stringsAsFactors=FALSE)
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase45_LR_expression_20261001'
d <- read.delim(gzfile(file.path(root,'tables/independent_LR_patient_P_QA_input.tsv.gz')))
out <- list()
for(ix in split(seq_len(nrow(d)),paste(d$context,d$pair_id,sep='|'))) {
 g<-d[ix,,drop=FALSE];v<-g$effect[is.finite(g$effect)];n<-length(v);p<-NA_real_
 if(n>=2 && !all(v==0) && n<=15) {
  signs<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));p<-mean(abs(as.vector(signs%*%v)/n)>=abs(mean(v))-1e-12)
 }
 expected<-g$expected_P[1];stopifnot(n==g$expected_N[1])
 if(n<=15)stopifnot((is.na(p)&&is.na(expected))||(!is.na(p)&&!is.na(expected)&&abs(p-expected)<1e-12))
 out[[length(out)+1]]<-data.frame(context=g$context[1],pair_id=g$pair_id[1],N=n,R_P=p,expected_P=expected,pass=TRUE)
}
write.table(do.call(rbind,out),file.path(root,'tables/independent_LR_patient_P_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
s<-read.delim(gzfile(file.path(root,'tables/independent_LR_site_score_reconstruction_input.tsv.gz')))
s$R_site_effect<-(log1p(s$ligand_target_CP10k)-log1p(s$ligand_reference_CP10k))*log1p(s$receptor_CP10k)
s$R_site_effect[!s$feature_available]<-NA_real_
lookup<-setNames(d$effect,paste(d$context,d$pair_id,d$donor,sep='|'));recon<-list()
for(ix in split(seq_len(nrow(s)),paste(s$context,s$pair_id,s$donor,sep='|'))) {
 g<-s[ix,,drop=FALSE];key<-paste(g$context[1],g$pair_id[1],g$donor[1],sep='|');v<-g$R_site_effect[is.finite(g$R_site_effect)];value<-if(length(v))mean(v)else NA_real_;expected<-unname(lookup[key])
 stopifnot((is.na(value)&&is.na(expected))||(!is.na(value)&&!is.na(expected)&&abs(value-expected)<1e-10))
 recon[[length(recon)+1]]<-data.frame(key=key,sites=nrow(g),R_effect=value,expected_effect=expected,pass=TRUE)
}
write.table(do.call(rbind,recon),file.path(root,'tables/independent_LR_score_reconstruction_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
f<-file.path(root,'tables/all_patient_LR_compatibility_statistics.tsv.gz');hdr<-names(read.delim(gzfile(f),nrows=0,check.names=FALSE));classes<-ifelse(hdr%in%c('P','FDR'),'numeric',ifelse(hdr=='family','character','NULL'));x<-read.delim(gzfile(f),colClasses=classes,check.names=FALSE);families<-list()
for(ix in split(seq_len(nrow(x)),x$family)) {
 g<-x[ix,,drop=FALSE];q<-p.adjust(ifelse(is.na(g$P),1,g$P),'BH');q[is.na(g$P)]<-NA_real_;stopifnot(identical(is.na(q),is.na(g$FDR)));err<-if(any(is.finite(q)))max(abs(q-g$FDR),na.rm=TRUE)else 0;stopifnot(err<1e-10)
 families[[length(families)+1]]<-data.frame(family=g$family[1],rows=nrow(g),valid_P=sum(!is.na(g$P)),max_error=err,pass=TRUE)
}
write.table(do.call(rbind,families),file.path(root,'tables/independent_LR_FDR_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(root,'INDEPENDENT_R_QA_COMPLETE'));cat('STAGE_D_LR_INDEPENDENT_QA_COMPLETE\n')