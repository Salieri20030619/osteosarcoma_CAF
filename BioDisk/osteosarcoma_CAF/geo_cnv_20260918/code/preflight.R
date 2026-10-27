.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages({library(infercnv);library(copykat);library(Matrix);library(data.table);library(jsonlite)})
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918'
a <- get('full.anno',envir=as.environment('package:copykat'))
s <- as.data.frame(fread(file.path(root,'review/input_summary.tsv')))
r <- list()
for (sm in s$sample[s$cnv_eligible]) {
 g <- as.data.frame(fread(file.path(root,'inputs',sm,'genes.tsv')))
 z <- a[a$ensembl_gene_id %in% g$ensembl_id,,drop=FALSE];z<-z[!duplicated(z$hgnc_symbol),,drop=FALSE]
 amb<-unique(z$ensembl_gene_id[duplicated(z$ensembl_gene_id)|duplicated(z$ensembl_gene_id,fromLast=TRUE)]);z<-z[!z$ensembl_gene_id %in% amb,,drop=FALSE]
 stopifnot(nrow(z)>5000)
 r[[sm]]<-data.frame(sample=sm,copykat_mapped_genes=nrow(z))
}
write.table(do.call(rbind,r),file.path(root,'review/copykat_mapping_preflight.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(R=as.character(getRversion()),infercnv=as.character(packageVersion('infercnv')),copykat=as.character(packageVersion('copykat')),samples=length(r)),file.path(root,'review/CNV_environment.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines('CNV_ENVIRONMENT_AND_MAPPING_PASS',file.path(root,'PREFLIGHT_PASS'))
cat('CNV_ENVIRONMENT_AND_MAPPING_PASS\n')
