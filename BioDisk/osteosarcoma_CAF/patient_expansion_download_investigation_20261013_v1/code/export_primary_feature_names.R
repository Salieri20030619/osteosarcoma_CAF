.libPaths(c("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib",.libPaths()))
suppressPackageStartupMessages(library(Seurat))
base <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1"
x <- qs2::qs_read(file.path(base,"downloads/OsteoCAR/patient_prim.qs2"),nthreads=2)
f <- data.frame(feature=rownames(x[["RNA"]]))
write.table(f,file.path(base,"tables/OsteoCAR/patient_prim_feature_names.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("GENE_MATCHES\n")
print(as.data.frame(f[grepl("RGS2|390957|PEF1|096854|35790|285280",f$feature),,drop=FALSE]))
cat("RNA_FEATURE_METADATA",paste(colnames(x[["RNA"]][[]]),collapse="|"),"\n")
