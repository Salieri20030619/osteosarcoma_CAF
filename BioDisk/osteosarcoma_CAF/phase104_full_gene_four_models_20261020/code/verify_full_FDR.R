P <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase104_full_gene_four_models_20261020"
models <- c("pooled_tumor","source_adjusted_tumor","pooled_normal_min20","source_adjusted_normal_min20")
results <- list()
for(m in models){
 x <- read.delim(gzfile(file.path(P,"tables",paste0(m,"_allgenes.tsv.gz"))),check.names=FALSE)
 valid <- is.finite(x$P)
 q <- p.adjust(x$P[valid],"BH")
 error <- max(abs(q-x$FDR[valid]))
 stopifnot(error<1e-12)
 results[[m]] <- data.frame(model=m,valid_tests=sum(valid),R_BH_max_error=error,FDR_significant=sum(q<0.05))
}
write.table(do.call(rbind,results),file.path(P,"resources/R_BH_VERIFICATION.tsv"),sep="\t",row.names=FALSE,quote=FALSE,fileEncoding="UTF-8")
print(as.data.frame(do.call(rbind,results)))
