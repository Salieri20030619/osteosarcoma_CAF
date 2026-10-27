suppressPackageStartupMessages(library(Seurat))
x<-Seurat::cc.genes.updated.2019
stopifnot(length(x$s.genes)>=40,length(x$g2m.genes)>=50)
write.table(data.frame(gene=c(x$s.genes,x$g2m.genes),phase=c(rep('S',length(x$s.genes)),rep('G2M',length(x$g2m.genes)))), 'remediation_20260917/resources/seurat_cc_updated2019.tsv',sep='\t',quote=FALSE,row.names=FALSE)
cat('SEURAT',as.character(packageVersion('Seurat')),'S',length(x$s.genes),'G2M',length(x$g2m.genes),'\n')
