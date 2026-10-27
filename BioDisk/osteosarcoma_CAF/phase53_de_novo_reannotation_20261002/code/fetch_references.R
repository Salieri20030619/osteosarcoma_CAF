options(stringsAsFactors=FALSE,timeout=600)
for(n in c("http_proxy","https_proxy","HTTP_PROXY","HTTPS_PROXY","ALL_PROXY","all_proxy"))Sys.unsetenv(n)
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SummarizedExperiment))
out <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002"
for(name in c("HPCA","BlueprintEncode")){
 message("FETCH_REFERENCE ",name)
 ref <- if(name=="HPCA") HumanPrimaryCellAtlasData() else BlueprintEncodeData()
 stopifnot(nrow(ref)>1000,ncol(ref)>10)
 saveRDS(ref,file.path(out,"resources",paste0(name,"reference.rds")))
 write.table(as.data.frame(colData(ref)),file.path(out,"resources",paste0(name,"metadata.tsv")),sep="\t",quote=FALSE,row.names=TRUE)
 write.table(data.frame(gene=rownames(ref)),file.path(out,"resources",paste0(name,"genes.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 cat(name,nrow(ref),ncol(ref),"\n")
 print(as.data.frame(table(label=ref$label.main)))
}
writeLines(capture.output(sessionInfo()),file.path(out,"resources","REFERENCE_R_SESSION.txt"))
writeLines("SUCCESS",file.path(out,"REFERENCES_COMPLETE"))
