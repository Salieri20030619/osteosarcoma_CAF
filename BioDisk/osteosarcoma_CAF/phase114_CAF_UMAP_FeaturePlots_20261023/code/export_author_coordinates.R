.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib',.libPaths()))
suppressPackageStartupMessages(library(Seurat))
p <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
o <- file.path(p,'phase114_CAF_UMAP_FeaturePlots_20261023')
x <- qs2::qs_read(file.path(p,'patient_expansion_download_investigation_20261013_v1/downloads/OsteoCAR/patient_mets.qs2'), nthreads=2)
cat('REDUCTIONS',names(x@reductions),'\n')
info <- list()
for(n in names(x@reductions)) {
 e <- Embeddings(x,n)
 info[[n]] <- list(dim=dim(e),assay=x@reductions[[n]]@assay.used)
 if(grepl('umap',n,ignore.case=TRUE)) {
 z <- data.frame(cell_id=paste0('NCI__',rownames(e)),e,check.names=FALSE)
 write.table(z,file.path(o,'tables',paste0('author_',n,'_coordinates.tsv')),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
 }
}
jsonlite::write_json(info,file.path(o,'resources/NCI_AUTHOR_REDUCTIONS.json'),pretty=TRUE,auto_unbox=TRUE)
