args<-commandArgs(TRUE);stopifnot(length(args)==2)
method<-args[1];sm<-args[2]
P<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF';O<-file.path(P,'phase91_lung_patient_validation_20261015');d<-file.path(O,'inputs',sm);out<-file.path(O,'cnv',method,sm)
.libPaths(c(file.path(P,'phase3_identity_20260914','Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages({library(Matrix);library(data.table);library(jsonlite)})
options(scipen=100);set.seed(20260922);dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'COMPLETE')))quit(status=0)
tryCatch({
 cells<-as.data.frame(fread(file.path(d,'cells.tsv')));genes<-as.data.frame(fread(file.path(d,'genes.tsv')))
 x<-as(readMM(gzfile(file.path(d,'counts.mtx.gz'))),'CsparseMatrix');dimnames(x)<-list(genes$ensembl,cells$tool_id)
 stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(abs(colSums(x)-cells$n_counts)<1e-6))
 params<-list(sample=sm,method=method,input_cells=ncol(x),input_genes=nrow(x),seed=20260922,project_numerical_entry_gates=FALSE,threads=4)
 if(method=='copykat'){
  suppressPackageStartupMessages(library(copykat));params$version<-as.character(packageVersion('copykat'));params$native_formals<-lapply(formals(copykat::copykat),function(z)paste(deparse(z),collapse=''))
  write_json(params,file.path(out,'parameters.json'),pretty=TRUE,auto_unbox=TRUE)
  setwd(out);ck<-copykat::copykat(rawmat=as.matrix(x),id.type='E',sam.name=sm,n.cores=4,output.seg=FALSE,plot.genes=FALSE)
  pr<-as.data.frame(ck$prediction);stopifnot(all(c('cell.names','copykat.pred') %in% names(pr)))
  saveRDS(ck,file.path(out,'copykat_result.rds'));calls<-as.character(pr$copykat.pred[match(cells$tool_id,pr$cell.names)]);calls[is.na(calls)]<-'not_returned_by_software'
  fwrite(data.table(cell_id=cells$cell_id,call=calls),file.path(out,'cell_results.tsv'),sep='\t')
 }else{
  suppressPackageStartupMessages(library(infercnv));params$version<-as.character(packageVersion('infercnv'));params$native_create_formals<-lapply(formals(infercnv::CreateInfercnvObject),function(z)paste(deparse(z),collapse=''));params$native_run_formals<-lapply(formals(infercnv::run),function(z)paste(deparse(z),collapse=''))
  refs<-if(any(cells$reference))'reference' else NULL
  params$reference_cells<-sum(cells$reference);params$baseline<-if(is.null(refs))'native_all_cell_mean_no_named_reference' else 'named_reference_from_source_annotations';params$recipe<-list(cutoff=.1,HMM=TRUE,HMM_type='i6',denoise=TRUE);params$cutoff_source<-'inferCNV official UMI/droplet recommendation';params$no_custom_signal_threshold<-TRUE
  write_json(params,file.path(out,'parameters.json'),pretty=TRUE,auto_unbox=TRUE)
  ic<-infercnv::CreateInfercnvObject(raw_counts_matrix=x,annotations_file=file.path(d,'annotations.tsv'),delim='\t',gene_order_file=file.path(P,'geo_cnv_20260918/resources/infercnv_gene_order_GRCh38_98.tsv'),ref_group_names=refs)
  ic<-infercnv::run(ic,cutoff=.1,out_dir=out,HMM=TRUE,HMM_type='i6',denoise=TRUE,num_threads=4,plot_steps=FALSE,save_rds=FALSE,save_final_rds=TRUE,write_expr_matrix=FALSE,output_format='png')
  gf<-list.files(out,pattern='^17_HMM.*cell_groupings$',full.names=TRUE);rf<-list.files(out,pattern='^HMM_CNV_predictions.*Pnorm_.*pred_cnv_regions.dat$',full.names=TRUE)
  calls<-rep('not_returned_by_software',nrow(cells))
  if(length(gf)==1L && length(rf)==1L){
   gr<-as.data.frame(fread(gf));regions<-as.data.frame(fread(rf));stopifnot(all(c('cell_group_name','cell') %in% names(gr)),all(c('cell_group_name','state') %in% names(regions)))
   # i6 state 3 is native neutral. No chromosome-count/amplitude/quantile rule is added.
   altered<-unique(regions$cell_group_name[regions$state!=3]);j<-match(cells$tool_id,gr$cell);ok<-!is.na(j);calls[ok]<-ifelse(gr$cell_group_name[j[ok]] %in% altered,'alteration_called','no_alteration_called')
  }
  fwrite(data.table(cell_id=cells$cell_id,call=calls),file.path(out,'cell_results.tsv'),sep='\t')
  if(!any(calls!='not_returned_by_software'))stop('Native HMM outputs not mapped; no normal calls assigned')
 }
 sink(file.path(out,'sessionInfo.txt'));print(sessionInfo());sink();writeLines('NATIVE_CNV_COMPLETE',file.path(out,'COMPLETE'));cat('NATIVE_CNV_COMPLETE',method,sm,'\n')
},error=function(e){writeLines(conditionMessage(e),file.path(out,'FAILED'));cat('NATIVE_CNV_FAILED',method,sm,conditionMessage(e),'\n');quit(status=1)})
