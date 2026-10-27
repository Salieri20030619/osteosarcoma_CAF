options(stringsAsFactors=FALSE)
root <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase50_F3_spatial_decomposition_20261002"
readtab <- function(n) read.delim(file.path(root,"tables",n), check.names=FALSE, na.strings=c("NA",""), quote="")
pat <- readtab("patient_effects.tsv")
st <- readtab("statistics.tsv")
p2 <- numeric(nrow(st));m2 <- numeric(nrow(st))
for(i in seq_len(nrow(st))){
 r <- st[i,]
 d <- pat$effect[pat$signature==r$signature & pat$model==r$model & pat$endpoint==r$endpoint]
 d <- d[is.finite(d)];stopifnot(length(d)==r$N)
 sums <- 0
 for(v in d) sums <- c(sums+v,sums-v)
 m2[i] <- mean(d)
 p2[i] <- mean(abs(sums/length(d))>=abs(mean(d))-1e-12)
}
q2 <- rep(NA_real_,nrow(st))
for(model in unique(st$model)){
 ii <- which(st$model==model & st$signature!="full_reference")
 stopifnot(length(ii)==12)
 q2[ii] <- p.adjust(p2[ii],method="BH")
}
statsQA <- data.frame(signature=st$signature,model=st$model,endpoint=st$endpoint,mean_error=abs(m2-st$mean_effect),P_error=abs(p2-st$P),FDR_error=abs(q2-st$FDR))
stopifnot(max(statsQA$mean_error)<1e-12,max(statsQA$P_error)<1e-12,max(statsQA$FDR_error,na.rm=TRUE)<1e-12)
write.table(statsQA,file.path(root,"tables/independent_R_statistics_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
sp <- read.delim(gzfile(file.path(root,"tables/spot_scores_for_QA.tsv.gz")),check.names=FALSE,na.strings=c("NA",""),quote="")
se <- readtab("sample_effects.tsv");se <- se[se$model=="full_adjusted",]
corrows <- vector("list",nrow(se))
for(i in seq_len(nrow(se))){
 r <- se[i,];d <- sp[sp$sample==r$sample,]
 x <- d[[r$signature]];y <- d[[r$endpoint]]
 covs <- as.matrix(d[,c("cancer","fibroblast","library_UMI","x","y")]);covs[,3] <- log1p(covs[,3])
 ok <- is.finite(x)&is.finite(y)&apply(is.finite(covs),1,all)
 rx <- rank(x[ok]);ry <- rank(y[ok]);cm <- cbind(1,apply(covs[ok,,drop=FALSE],2,rank))
 ex <- lm.fit(cm,rx)$residuals;ey <- lm.fit(cm,ry)$residuals
 value <- cor(ex,ey)
 corrows[[i]] <- data.frame(sample=r$sample,signature=r$signature,endpoint=r$endpoint,R_effect=value,Python_effect=r$effect,abs_error=abs(value-r$effect))
}
cq <- do.call(rbind,corrows)
write.table(cq,file.path(root,"tables/independent_R_correlations_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(nrow(cq)==364,max(cq$abs_error)<1e-10)
out <- data.frame(check=c("patient_statistic_rows","primary_sample_correlations","max_mean_error","max_P_error","max_FDR_error","max_correlation_error"),value=c(nrow(st),nrow(cq),max(statsQA$mean_error),max(statsQA$P_error),max(statsQA$FDR_error,na.rm=TRUE),max(cq$abs_error)))
write.table(out,file.path(root,"delivery/INDEPENDENT_R_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines("SUCCESS",file.path(root,"INDEPENDENT_R_QA_COMPLETE"))
print(as.data.frame(out))
