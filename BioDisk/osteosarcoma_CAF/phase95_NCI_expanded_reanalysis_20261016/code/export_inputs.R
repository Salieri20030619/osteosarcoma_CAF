.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib',.libPaths()))
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(jsonlite)})
P<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
O<-file.path(P,'phase95_NCI_expanded_reanalysis_20261016')
src<-file.path(P,'patient_expansion_download_investigation_20261013_v1/downloads/OsteoCAR/patient_mets.qs2')
cat('READ_SOURCE\n');flush.console()
obj<-qs2::qs_read(src,nthreads=2)
md<-obj@meta.data;sel<-md$data_source=='NCI_POB';stopifnot(sum(sel)==123255,length(unique(md$sample_name[sel]))==19)
x<-SeuratObject::GetAssayData(obj,assay='RNA',layer='counts')
stopifnot(ncol(x)==nrow(md),identical(colnames(x),rownames(md)),!anyDuplicated(rownames(x)),all(x@x>=0),all(x@x==floor(x@x)))
manifest<-list()
for(sm in sort(unique(as.character(md$sample_name[sel])))){
 ix<-which(sel&md$sample_name==sm);y<-x[,ix,drop=FALSE];m<-md[ix,,drop=FALSE]
 stopifnot(all(m$location=='lung'),all(m$method=='single_nucleus'),all(m$organism=='human'))
 name<-paste0('NCI__',sm);d<-file.path(O,'inputs',name);dir.create(d,recursive=TRUE,showWarnings=FALSE)
 m$author_cell_id<-rownames(m);m$cell_id<-paste0('NCI__',rownames(m));m$sample<-name
 m$tool_id<-sprintf('%s_c%06d',name,seq_len(nrow(m))-1);m$n_counts<-as.numeric(colSums(y));m$n_genes<-as.numeric(colSums(y>0));m$site<-'Lung'
 write.table(m,file.path(d,'author_cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')
 write.table(data.frame(feature_id=rownames(y),symbol=rownames(y)),file.path(d,'genes.tsv'),sep='\t',quote=FALSE,row.names=FALSE,fileEncoding='UTF-8')
 raw<-file.path(d,'counts.mtx');writeMM(y,raw)
 stopifnot(system2('gzip',c('-n','-f',shQuote(raw)))==0)
 check<-as(readMM(gzfile(paste0(raw,'.gz'))),'CsparseMatrix');stopifnot(identical(dim(check),dim(y)),all(check==y))
 manifest[[sm]]<-list(sample=name,author_sample=sm,cells=ncol(y),features=nrow(y),source='OsteoCAR patient_mets RNA counts; processed available-cell object',author_CAF=sum(m$Ann_Level0=='Host'&m$Ann_Level2=='CAFs'),roundtrip=TRUE)
 writeLines('EXPORTED_AVAILABLE_INPUT',file.path(d,'EXPORT_COMPLETE'));cat('EXPORTED',name,ncol(y),'\n');flush.console()
 rm(y,m,check);gc(FALSE)
}
write_json(manifest,file.path(O,'resources/EXPORTED_INPUTS.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines('EXPORT_COMPLETE',file.path(O,'resources/EXPORT_COMPLETE'))
