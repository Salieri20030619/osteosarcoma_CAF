# SYNTAX_CHECK: Rscript --vanilla check_stageC_R.R
# RUN_COMMAND: Rscript --vanilla stageC_independent_QA.R
# EXPECTED_SUCCESS_MARKER: STAGE_C1_INDEPENDENT_R_QA_COMPLETE
options(stringsAsFactors=FALSE)
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase43_full_transcriptome_discovery_20261001'
d <- read.delim(gzfile(file.path(root,'tables/independent_R_patient_effect_QA_input.tsv.gz')),check.names=FALSE)
key <- paste(d$context,d$feature,sep='|')
out <- list()
for (ix in split(seq_len(nrow(d)),key)) {
 g <- d[ix,,drop=FALSE];v <- g$effect[is.finite(g$effect)];n <- length(v)
 if (n<2 || all(v==0)) p <- NA_real_ else if(n<=15) {
  signs <- as.matrix(expand.grid(rep(list(c(-1,1)),n)))
  null <- abs(as.vector(signs %*% v)/n)
  p <- mean(null >= abs(mean(v))-1e-12)
 } else p <- NA_real_
 expected <- g$expected_P[1]
 if(n<=15) stopifnot((is.na(p)&&is.na(expected))||(!is.na(p)&&!is.na(expected)&&abs(p-expected)<1e-12))
 stopifnot(n==g$expected_N[1])
 if(n>0) stopifnot(abs(mean(v)-g$expected_mean[1])<1e-10,abs(median(v)-g$expected_median[1])<1e-10)
 out[[length(out)+1]] <- data.frame(context=g$context[1],feature=g$feature[1],N=n,Python_P=expected,R_P=p,exact_checked=n<=15,pass=TRUE)
}
write.table(do.call(rbind,out),file.path(root,'tables/independent_R_patient_statistics_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
fdrout <- list()
files <- list.files(file.path(root,'tables'),pattern='__full_gene_statistics.tsv.gz$',full.names=TRUE)
for (f in files) {
 header <- names(read.delim(gzfile(f),nrows=0,check.names=FALSE))
 classes <- ifelse(header %in% c('P','FDR'),'numeric',ifelse(header=='family','character','NULL'))
 x <- read.delim(gzfile(f),colClasses=classes,check.names=FALSE)
 for(ix in split(seq_len(nrow(x)),x$family)) {
  g <- x[ix,,drop=FALSE];q <- p.adjust(ifelse(is.na(g$P),1,g$P),method='BH');q[is.na(g$P)] <- NA_real_
  stopifnot(identical(is.na(q),is.na(g$FDR)))
  err <- if(any(is.finite(q))) max(abs(q-g$FDR),na.rm=TRUE) else 0
  stopifnot(err<1e-10)
  fdrout[[length(fdrout)+1]] <- data.frame(file=basename(f),family=g$family[1],rows=nrow(g),valid_P=sum(!is.na(g$P)),max_abs_error=err,pass=TRUE)
 }
 rm(x);gc(verbose=FALSE)
}
write.table(do.call(rbind,fdrout),file.path(root,'tables/independent_R_FDR_family_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(root,'INDEPENDENT_R_QA_COMPLETE'))
cat('STAGE_C1_INDEPENDENT_R_QA_COMPLETE\n')