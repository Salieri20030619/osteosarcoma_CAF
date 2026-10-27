args<-commandArgs(TRUE);stopifnot(length(args)==2L)
kind<-args[1];sm<-args[2]
P<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF';B<-file.path(P,'phase26_client_reanalysis_20260922');O<-file.path(P,'phase31_customer_completion_20260926')
.libPaths(c(file.path(P,'phase3_identity_20260914/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(Matrix);library(SingleCellExperiment);library(jsonlite)})
d<-file.path(O,kind,sm);dir.create(d,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(d,'COMPLETE')))stop('Existing completed output; controller must reuse it')
tryCatch({
 cells<-read.delim(file.path(B,'inputs',sm,'cells.tsv'),check.names=FALSE)
 genes<-read.delim(file.path(B,'inputs',sm,'genes.tsv'),check.names=FALSE)
 x<-as(readMM(gzfile(file.path(B,'inputs',sm,'counts.mtx.gz'))),'CsparseMatrix')
 stopifnot(ncol(x)==nrow(cells),nrow(x)==nrow(genes),all(x@x>=0),all(x@x==floor(x@x)))
 dimnames(x)<-list(genes$ensembl,cells$cell_id)
 sce<-SingleCellExperiment(list(counts=x));set.seed(20260926)
 params<-list(sample=sm,method=kind,input_cells=ncol(x),input_genes=nrow(x),project_cell_filter=FALSE,seed=20260926,software_defaults_except_seed_and_computational_workers=TRUE)
 if(kind=='decontx'){
  suppressPackageStartupMessages(library(decontX))
  params$version<-as.character(packageVersion('decontX'));params$reference<-'no empty-droplet background; native automatic broad clustering; no invented contamination acceptance cutoff'
  write_json(params,file.path(d,'parameters.json'),auto_unbox=TRUE,pretty=TRUE)
  sce<-decontX::decontX(sce,seed=20260926)
  y<-assay(sce,'decontXcounts');cont<-colData(sce)$decontX_contamination
  stopifnot(identical(dim(y),dim(x)),identical(colnames(y),cells$cell_id),all(is.finite(cont)),all(cont>=0 & cont<=1),all(is.finite(y@x)),all(y@x>=0))
  write.table(data.frame(cell_id=cells$cell_id,contamination=cont,original_UMI=colSums(x),corrected_UMI=colSums(y)),file.path(d,'cell_contamination.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeMM(y,file.path(d,'corrected_counts.mtx'))
  stopifnot(system2('gzip',c('-f',shQuote(file.path(d,'corrected_counts.mtx'))))==0)
 } else if(kind=='doublet'){
  suppressPackageStartupMessages({library(scDblFinder);library(BiocParallel)})
  params$version<-as.character(packageVersion('scDblFinder'));params$project_doublet_rate_override<-FALSE
  write_json(params,file.path(d,'parameters.json'),auto_unbox=TRUE,pretty=TRUE)
  sce<-scDblFinder::scDblFinder(sce,BPPARAM=BiocParallel::SerialParam())
  stopifnot(identical(colnames(sce),cells$cell_id))
  write.table(data.frame(cell_id=cells$cell_id,scDblFinder_score=colData(sce)$scDblFinder.score,scDblFinder_class=colData(sce)$scDblFinder.class),file.path(d,'doublet_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 } else stop('unknown method')
 sink(file.path(d,'sessionInfo.txt'));print(sessionInfo());sink()
 writeLines('SOFTWARE_OUTPUT_AND_CELL_IDS_VERIFIED',file.path(d,'COMPLETE'))
 cat('COMPLETE',kind,sm,'\n')
},error=function(e){writeLines(conditionMessage(e),file.path(d,'FAILED'));cat('FAILED',kind,sm,conditionMessage(e),'\n');quit(status=1)})
