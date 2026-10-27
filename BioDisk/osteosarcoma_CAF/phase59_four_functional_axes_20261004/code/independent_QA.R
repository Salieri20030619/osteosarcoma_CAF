options(stringsAsFactors=FALSE)
o <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase59_four_functional_axes_20261004'
readtab <- function(p) read.delim(p, check.names=FALSE, na.strings=c('','NA','nan'))
s <- readtab(file.path(o,'tables/axis_statistics.tsv'))
e <- readtab(gzfile(file.path(o,'tables/patient_axis_effects.tsv.gz')))
g <- readtab(gzfile(file.path(o,'tables/marker_statistics.tsv.gz')))
ge <- readtab(gzfile(file.path(o,'tables/patient_marker_effects.tsv.gz')))
key <- function(d,item) do.call(paste,c(d[c('dataset','matrix','context','contrast',item)],sep='|'))
check <- function(s,e,item) {
 sk <- key(s,item); ek <- key(e,item); stopifnot(!anyDuplicated(sk)); stopifnot(!anyDuplicated(paste(ek,e$donor_record)))
 z <- split(e,ek); vals <- vector('list',nrow(s)); errors <- rep(NA_real_,nrow(s)); counts <- integer(nrow(s))
 for(i in seq_len(nrow(s))) {
  d <- z[[sk[i]]]; v <- if(is.null(d)) numeric(0) else d$effect[order(d$donor_record)]; v <- v[is.finite(v)];vals[[i]] <- v; counts[i]<-length(v)
  stopifnot(length(v)==s$N[i]); if(length(v))stopifnot(abs(mean(v)-s$mean_delta[i])<1e-9)
  stopifnot(sum(v>0)==s$positive_patients[i],sum(v<0)==s$negative_patients[i])
 }
 for(n in sort(unique(counts[counts>=2]))) {
  ii <- which(counts==n & !is.na(s$P)); if(!length(ii))next
  if(n<=18) {
   integers <- 0:(2^n-1); signs <- sapply(0:(n-1),function(k) ifelse(bitwAnd(bitwShiftR(integers,k),1L)==1L,1,-1)); exact <- TRUE
  } else { signs <- as.matrix(read.table(gzfile(file.path(o,'QA',paste0('MC_signs_N',n,'.tsv.gz'))),header=FALSE));exact <- FALSE }
  for(st in seq(1,length(ii),by=32)) {
   jj <- ii[st:min(st+31,length(ii))]; mat <- do.call(cbind,vals[jj]); obs <- abs(colMeans(mat)); null <- abs(signs %*% mat/n); hits <- colSums(sweep(null,2,obs-1e-12,'>=')); pv <- if(exact) hits/nrow(signs) else (hits+1)/(nrow(signs)+1);errors[jj] <- abs(pv-s$P[jj])
  }
 }
 for(i in which(is.na(s$P))) stopifnot(length(vals[[i]])<2 || all(vals[[i]]==0))
 stopifnot(max(errors,na.rm=TRUE)<1e-9)
 ferr <- numeric(0)
 for(f in unique(s$family)) {
  ii <- which(s$family==f); pv <- s$P[ii]; pv[is.na(pv)]<-1; q<-p.adjust(pv,'BH');keep<-!is.na(s$P[ii]);ferr<-c(ferr,abs(q[keep]-s$FDR[ii][keep]));stopifnot(all(is.na(s$FDR[ii][!keep])))
 }
 stopifnot(max(ferr,na.rm=TRUE)<1e-9)
 data.frame(type=item,rows=nrow(s),valid_P=sum(!is.na(s$P)),families=length(unique(s$family)),max_P_error=max(errors,na.rm=TRUE),max_FDR_error=max(ferr,na.rm=TRUE))
}
result <- rbind(check(s,e,'program'),check(g,ge,'gene'));write.table(result,file.path(o,'QA/independent_R_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(as.data.frame(result));writeLines('SUCCESS',file.path(o,'QA/INDEPENDENT_R_COMPLETE'))
