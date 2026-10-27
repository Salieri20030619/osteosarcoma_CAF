# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' run_scevan.R
# RUN_COMMAND: Rscript --vanilla run_scevan.R <sample_short>
# EXPECTED_SUCCESS_MARKER: SCEVAN_SAMPLE_COMPLETE
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1);s<-args[1]
root<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917'
.libPaths(c('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/SCEVAN_Rlib','/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/Rlib',.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages({library(SCEVAN);library(Matrix);library(data.table);library(jsonlite)})
out<-file.path(root,'scevan',s);dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'COMPLETE'))) {cat('SCEVAN_CACHED\n');quit(status=0)}
status<-function(state,extra=list()) write_json(c(list(method='SCEVAN',sample=s,state=state,updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),extra),file.path(out,'status.json'),auto_unbox=TRUE,pretty=TRUE)
status('running');options(scipen=100)
tryCatch({
  input<-file.path(root,'inputs',s)
  cells<-as.data.frame(fread(file.path(input,'cells.tsv')));genes<-as.data.frame(fread(file.path(input,'genes.tsv')))
  mat<-as(readMM(gzfile(file.path(input,'counts.mtx.gz'))),'CsparseMatrix');dimnames(mat)<-list(genes$ensembl_id,cells$tool_id)
  stopifnot(all(mat@x>=0),all(mat@x==floor(mat@x)),all(abs(colSums(mat)-cells$n_counts)<1e-6))
  set.seed(270231)
  # Use the package's actual preprocessing/classification stages. No fibroblast normal signature,
  # no forced FIXED_NORMAL_CELLS, and no unnecessary subclone analysis after classification.
  proc<-getFromNamespace('preprocessingMtx','SCEVAN')(as.matrix(mat),sample=s,ngenes_chr=5,perc_genes=.1,par_cores=4,findConfident=FALSE,SCEVANsignatures=FALSE,organism='human',output_dir=out)
  available<-colnames(proc$count_mtx_norm)
  ri<-intersect(cells$tool_id[cells$cnv_role=='ref_immune'],available)
  re<-intersect(cells$tool_id[cells$cnv_role=='ref_endo'],available)
  set.seed(270231)
  if(length(ri)>=20 && length(re)>=20) {
    n<-min(length(ri),length(re));refs<-c(sample(ri,n),sample(re,n));strategy<-'same-sample balanced immune/endo'
  } else {
    refs<-if(length(ri)>=30)ri else if(length(re)>=30)re else character(0)
    stopifnot(length(refs)>=30);n<-length(refs);strategy<-'same-sample single adequate normal lineage; unavailable second lineage not invented'
  }
  cells$cnv_role[grepl('^ref_',cells$cnv_role)&!cells$tool_id %in% refs]<-'check_ref'
  fwrite(data.table(tool_id=refs),file.path(out,'reference_cells.tsv'),sep='\t')
  res<-getFromNamespace('classifyTumorCells','SCEVAN')(proc$count_mtx_norm,proc$count_mtx_annot,sample=s,par_cores=4,norm_cell_names=refs,SEGMENTATION_CLASS=TRUE,SMOOTH=TRUE,beta_vega=.5,FIXED_NORMAL_CELLS=FALSE,output_dir=out)
  saveRDS(res,file.path(out,'classification_result.rds'))
  fwrite(as.data.table(proc$count_mtx_annot),file.path(out,'gene_positions_used.tsv'),sep='\t')
  cna<-as.data.frame(res$CNAmat,check.names=FALSE)
  ids<-intersect(cells$tool_id,colnames(cna));stopifnot(length(ids)>=15,nrow(cna)>=2000)
  values<-as.matrix(cna[,ids,drop=FALSE]);stopifnot(all(is.finite(values)))
  results<-cells[,c('tool_id','cell_id','sample','short','patient_id','barcode','cluster','phase2_label','lineage_phase3','cnv_role','n_counts','n_genes')]
  results$method<-'SCEVAN';results$call<-'not_evaluable';results$signal_strength<-NA_real_
  results$call[results$tool_id %in% ids]<-'normal';results$call[results$tool_id %in% res$tum_cells]<-'tumor'
  results$signal_strength[match(ids,results$tool_id)]<-colMeans(abs(values))
  results$reference_used<-results$tool_id %in% refs
  fwrite(results,file.path(out,'cell_results.tsv'),sep='\t')
  bins<-data.frame(chr=as.character(cna[[2]]),position=as.numeric(cna[[3]]));stopifnot(all(bins$chr %in% as.character(1:22)))
  saveRDS(list(bins=bins,matrix=values),file.path(out,'window_profiles.rds'))
  params<-list(method='SCEVAN',version=as.character(packageVersion('SCEVAN')),genome='GRCh38',annotation='EnsDb.Hsapiens.v86 internal table',input_cells=ncol(mat),evaluable_cells=length(ids),genes=nrow(cna),ngenes_chr=5,perc_genes=.1,beta_vega=.5,reference_per_lineage=n,reference=strategy,automatic_normal_signature_search=FALSE,FIXED_NORMAL_CELLS=FALSE,preprocessing_internal_seed=123,classification_internal_seed=1,reference_sampling_seed=270231,par_cores=4,stages='official preprocessingMtx + classifyTumorCells; classification-only')
  write_json(params,file.path(out,'parameters.json'),auto_unbox=TRUE,pretty=TRUE)
  sink(file.path(out,'sessionInfo.txt'));print(sessionInfo());sink()
  status('complete',list(calls=as.list(table(results$call))))
  writeLines(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),file.path(out,'COMPLETE'));cat('SCEVAN_SAMPLE_COMPLETE',s,'\n')
},error=function(e){status('failed',list(error=conditionMessage(e)));writeLines(conditionMessage(e),file.path(out,'FAILED'));cat('SCEVAN_SAMPLE_FAILED',conditionMessage(e),'\n');quit(status=1)})
