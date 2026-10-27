suppressPackageStartupMessages(library(Seurat))
p <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase97_normal_background_20261020"
x <- readRDS(file.path(p,"inputs/bone_msc_verified.rds"))
writeLines(capture.output(str(x,max.level=1)),file.path(p,"resources/bone_object_structure.txt"))
m <- x[[]]
write.table(data.frame(cell_id=rownames(m),m,check.names=FALSE),file.path(p,"tables/bone_metadata.tsv"),sep="\t",fileEncoding="UTF-8",quote=FALSE,row.names=FALSE)
cat("ASSAYS",names(x@assays),"\n")
print(colnames(m))
for (k in names(m)) if(length(unique(m[[k]]))<60) { cat("FIELD",k,"\n"); print(table(m[[k]],useNA="ifany")) }
cat("RNA features",nrow(x[["RNA"]]),"\n")
