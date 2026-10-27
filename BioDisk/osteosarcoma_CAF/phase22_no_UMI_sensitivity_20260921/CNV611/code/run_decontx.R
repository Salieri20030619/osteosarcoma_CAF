args<-commandArgs(TRUE);stopifnot(length(args)==1);sample<-args[1]
root<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase22_no_UMI_sensitivity_20260921/CNV611'
suppressPackageStartupMessages({library(decontX);library(Matrix);library(SingleCellExperiment);library(jsonlite)})
input<-file.path(root,'inputs',sample);out<-file.path(root,'decontx',sample);dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(input,'COMPLETE')))
cells<-read.delim(file.path(input,'cells.tsv'));genes<-read.delim(file.path(input,'genes.tsv'))
mat<-as(readMM(gzfile(file.path(input,'counts.mtx.gz'))),'CsparseMatrix');dimnames(mat)<-list(genes$gene_id,cells$cell_id)
stopifnot(ncol(mat)==nrow(cells),nrow(mat)==nrow(genes),all(mat@x>=0),all(mat@x==floor(mat@x)))
sce<-SingleCellExperiment(list(counts=mat));z<-as.integer(factor(cells$ambient_group));stopifnot(length(unique(z))>=2)
sce<-decontX::decontX(sce,z=z,seed=20260919,verbose=TRUE)
corrected<-assay(sce,'decontXcounts');cont<-colData(sce)$decontX_contamination
stopifnot(all(is.finite(cont)),all(cont>=0 & cont<=1),identical(dim(corrected),dim(mat)),all(is.finite(corrected@x)),all(corrected@x>=0))
write.table(data.frame(cell_id=cells$cell_id,cluster=cells$cluster,contamination=cont,original_UMI=colSums(mat),corrected_UMI=colSums(corrected)),file.path(out,'cell_contamination.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
Matrix::writeMM(corrected,file.path(out,'corrected_counts.mtx'));stopifnot(system2('gzip',c('-f',shQuote(file.path(out,'corrected_counts.mtx'))))==0)
saveRDS(metadata(sce)$decontX,file.path(out,'model_metadata.rds'))
write_json(list(sample=sample,cells=ncol(mat),genes=nrow(mat),package_version=as.character(packageVersion('decontX')),background='none; same-sample frozen broad lineage groups; sensitivity correction',median_contamination=median(cont),seed=20260919),file.path(out,'parameters.json'),auto_unbox=TRUE,pretty=TRUE)
sink(file.path(out,'sessionInfo.txt'));print(sessionInfo());sink()
writeLines('DECONTX_SAMPLE_COMPLETE',file.path(out,'COMPLETE'));cat('DECONTX_SAMPLE_COMPLETE',sample,'\n')
