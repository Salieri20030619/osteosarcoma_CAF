# SYNTAX_CHECK: Rscript --vanilla check_stageD_R.R
# RUN_COMMAND: Rscript --vanilla stageD_independent_QA.R
# EXPECTED_SUCCESS_MARKER: STAGE_D_INDEPENDENT_R_QA_COMPLETE
options(stringsAsFactors=FALSE)
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase44_patient_mechanism_20261001'
flip <- function(v) {
 v <- v[is.finite(v)];n <- length(v)
 if(n<2 || all(v==0))return(NA_real_)
 if(n>15)return(NA_real_)
 signs <- as.matrix(expand.grid(rep(list(c(-1,1)),n)))
 mean(abs(as.vector(signs %*% v)/n)>=abs(mean(v))-1e-12)
}
check_equal <- function(x,y,tol=1e-10) {
 if(is.na(x)&&is.na(y))return(TRUE)
 stopifnot(is.finite(x),is.finite(y),abs(x-y)<tol);TRUE
}
tf <- read.delim(gzfile(file.path(root,'tables/independent_TF_QA_patient_scores.tsv.gz')))
out <- list()
for(ix in split(seq_len(nrow(tf)),paste(tf$context,tf$regulator,sep='|'))) {
 g <- tf[ix,,drop=FALSE];v <- g$score[is.finite(g$score)];n <- length(v);p <- flip(v)
 stopifnot(n==g$expected_N[1]);if(n>0)check_equal(mean(v),g$expected_mean[1]);if(n<=15)check_equal(p,g$expected_P[1],1e-12)
 out[[length(out)+1]] <- data.frame(context=g$context[1],regulator=g$regulator[1],N=n,R_P=p,expected_P=g$expected_P[1],exact_checked=n<=15,pass=TRUE)
}
write.table(do.call(rbind,out),file.path(root,'tables/independent_TF_patient_P_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
target <- read.delim(gzfile(file.path(root,'tables/independent_TF_target_reconstruction_input.tsv.gz')))
score_lookup <- setNames(tf$score,paste(tf$context,tf$regulator,tf$donor,sep='|'));recon <- list()
for(ix in split(seq_len(nrow(target)),paste(target$context,target$regulator,target$donor,sep='|'))) {
 g <- target[ix,,drop=FALSE];ok <- is.finite(g$gene_effect);den <- sum(abs(g$weight[ok]));value <- if(den>0)sum(g$weight[ok]*g$gene_effect[ok])/den else NA_real_
 key <- paste(g$context[1],g$regulator[1],g$donor[1],sep='|');expected <- unname(score_lookup[key]);check_equal(value,expected)
 recon[[length(recon)+1]] <- data.frame(key=key,targets=nrow(g),measured_effects=sum(ok),R_score=value,expected_score=expected,pass=TRUE)
}
write.table(do.call(rbind,recon),file.path(root,'tables/independent_TF_signed_score_reconstruction_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
cl <- read.delim(file.path(root,'tables/independent_clinical_paired_QA_input.tsv'))
paired <- list()
if(nrow(cl)>0)for(ix in split(seq_len(nrow(cl)),cl$comparison_key)) {
 g <- cl[ix,,drop=FALSE];n <- length(g$paired_effect);p <- flip(g$paired_effect);stopifnot(n==g$expected_N[1]);if(n<=15)check_equal(p,g$expected_P[1],1e-12)
 paired[[length(paired)+1]] <- data.frame(key=g$comparison_key[1],N=n,R_P=p,expected_P=g$expected_P[1],exact_checked=n<=15,pass=TRUE)
}
if(length(paired))write.table(do.call(rbind,paired),file.path(root,'tables/independent_clinical_paired_P_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
families <- list()
for(name in c('all_TF_signed_patient_statistics.tsv.gz','all_clinical_fraction_statistics.tsv','adjusted_clinical_fraction_statistics.tsv')) {
 f <- file.path(root,'tables',name);conn <- if(grepl('gz$',name))gzfile(f,open='rt') else file(f,open='rt');hdr <- names(read.delim(conn,nrows=0,check.names=FALSE));close(conn)
 classes <- ifelse(hdr %in% c('P','FDR'),'numeric',ifelse(hdr=='family','character','NULL'));conn <- if(grepl('gz$',name))gzfile(f,open='rt') else file(f,open='rt');d <- read.delim(conn,colClasses=classes,check.names=FALSE);close(conn)
 for(ix in split(seq_len(nrow(d)),d$family)) {
  g <- d[ix,,drop=FALSE];q <- p.adjust(ifelse(is.na(g$P),1,g$P),'BH');q[is.na(g$P)] <- NA_real_;stopifnot(identical(is.na(q),is.na(g$FDR)));err <- if(any(is.finite(q)))max(abs(q-g$FDR),na.rm=TRUE)else 0;stopifnot(err<1e-10)
  families[[length(families)+1]] <- data.frame(file=name,family=g$family[1],rows=nrow(g),valid_P=sum(!is.na(g$P)),max_error=err,pass=TRUE)
 }
}
write.table(do.call(rbind,families),file.path(root,'tables/independent_stageD_FDR_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('SUCCESS',file.path(root,'INDEPENDENT_R_QA_COMPLETE'))
cat('STAGE_D_INDEPENDENT_R_QA_COMPLETE\n')