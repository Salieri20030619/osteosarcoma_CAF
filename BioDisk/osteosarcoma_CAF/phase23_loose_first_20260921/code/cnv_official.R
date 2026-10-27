args<-commandArgs(TRUE);stopifnot(length(args)==2)
method<-args[1];sm<-args[2]
P<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
O<-file.path(P,'phase23_loose_first_20260921','CNV_official_inputs')
.libPaths(c(file.path(P,'phase3_identity_20260914','Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages({library(Matrix);library(data.table);library(jsonlite)})
options(scipen=100);set.seed(20260921)
d<-file.path(O,'inputs',sm);out<-file.path(O,'results',method,sm);dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'COMPLETE')))quit(status=0)
tryCatch({
 cells<-as.data.frame(fread(file.path(d,'cells.tsv')));genes<-as.data.frame(fread(file.path(d,'genes.tsv')))
 x<-as(readMM(gzfile(file.path(d,'counts.mtx.gz'))),'CsparseMatrix');dimnames(x)<-list(genes$ensembl,cells$tool_id)
 stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(abs(colSums(x)-cells$n_counts)<1e-6))
 params<-list(method=method,sample=sm,input_cells=ncol(x),input_genes=nrow(x),no_project_entry_thresholds=TRUE,provider_processed_input=TRUE,seed=20260921)
 if(method=='copykat'){
  suppressPackageStartupMessages(library(copykat));params$package_version<-as.character(packageVersion('copykat'));params$native_formals<-lapply(formals(copykat::copykat),function(z)paste(deparse(z),collapse=''))
  # Leave native thresholds and normal-baseline inference at package defaults.
  setwd(out);ck<-copykat::copykat(rawmat=as.matrix(x),id.type='E',sam.name=sm,n.cores=4,output.seg=FALSE,plot.genes=FALSE)
  saveRDS(ck,file.path(out,'copykat_result.rds'));pr<-as.data.frame(ck$prediction);stopifnot(all(c('cell.names','copykat.pred') %in% names(pr)))
  calls<-as.character(pr$copykat.pred[match(cells$tool_id,pr$cell.names)]);calls[is.na(calls)]<-'not_returned_by_software'
  fwrite(data.table(cell_id=cells$cell_id,call=calls),file.path(out,'cell_results.tsv'),sep='\t')
 }else{
  suppressPackageStartupMessages(library(infercnv));params$package_version<-as.character(packageVersion('infercnv'));params$create_formals<-lapply(formals(infercnv::CreateInfercnvObject),function(z)paste(deparse(z),collapse=''));params$run_formals<-lapply(formals(infercnv::run),function(z)paste(deparse(z),collapse=''))
  refs<-if(any(cells$reference))'reference' else NULL;params$reference_cells<-sum(cells$reference);params$reference_source<-'provider is_infercnv_reference';params$HMM<-TRUE;params$HMM_type<-'i6';params$cutoff<-.1;params$denoise<-TRUE
  ic<-infercnv::CreateInfercnvObject(raw_counts_matrix=x,annotations_file=file.path(d,'annotations.tsv'),delim='\t',gene_order_file=file.path(P,'geo_cnv_20260918/resources/infercnv_gene_order_GRCh38_98.tsv'),ref_group_names=refs)
  # Official documented 10x recipe; all remaining filtering/HMM thresholds use native defaults.
  ic<-infercnv::run(ic,cutoff=.1,out_dir=out,HMM=TRUE,HMM_type='i6',denoise=TRUE,num_threads=4,output_format='png',plot_steps=FALSE,save_rds=FALSE,save_final_rds=TRUE,write_expr_matrix=FALSE)
  saveRDS(ic,file.path(out,'infercnv_result.rds'))
  # Continuous descriptive quantity only: no custom low/high threshold.
  score<-colMeans(abs(log2(pmax(ic@expr.data,1e-3))))
  fwrite(data.table(cell_id=cells$cell_id,mean_abs_log2=unname(score[match(cells$tool_id,names(score))])),file.path(out,'cell_continuous_signal.tsv'),sep='\t')
 }
 write_json(params,file.path(out,'parameters.json'),auto_unbox=TRUE,pretty=TRUE);sink(file.path(out,'sessionInfo.txt'));print(sessionInfo());sink()
 writeLines('NATIVE_METHOD_COMPLETE',file.path(out,'COMPLETE'));cat('NATIVE_METHOD_COMPLETE',method,sm,'\n')
},error=function(e){writeLines(conditionMessage(e),file.path(out,'FAILED'));cat('NATIVE_METHOD_FAILED',method,sm,conditionMessage(e),'\n');quit(status=1)})
