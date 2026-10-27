# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' export_cnv_plot_data.R
# RUN_COMMAND: Rscript --vanilla export_cnv_plot_data.R
# EXPECTED_SUCCESS_MARKER: CNV_PLOT_DATA_EXPORTED
suppressPackageStartupMessages(library(data.table))
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914'
stopifnot(file.exists(file.path(root,'identity','COMPLETE')))
out<-file.path(root,'identity','plot_data');dir.create(out,showWarnings=FALSE)
set.seed(270231)
for(method in c('infercnv','copykat','copykat_balanced','scevan')) for(s in c('S0058','S0059','S0217','S0218','SC069','SC072','SC073')) {
  if(!file.exists(file.path(root,method,s,'COMPLETE'))) next
  obj<-readRDS(file.path(root,method,s,'window_profiles.rds'))
  cells<-as.data.frame(fread(file.path(root,method,s,'cell_results.tsv')))
  splitids<-split(cells$tool_id,interaction(cells$cnv_role,cells$cluster,drop=TRUE))
  ids<-unlist(lapply(splitids,function(x) if(length(x)>25) sample(x,25) else x),use.names=FALSE)
  ids<-intersect(ids,colnames(obj$matrix));stopifnot(length(ids)>0)
  field<-intersect(c('chr','chrom','chromosome_name','chromosome'),names(obj$bins))[1]
  stopifnot(!is.na(field))
  chr<-sub('^chr','',as.character(obj$bins[[field]]));v<-obj$matrix[,ids,drop=FALSE]
  ans<-list();meta<-list();n<-0L
  for(cc in as.character(1:22)) {
    ix<-which(chr==cc)
    if(!length(ix)) next
    groups<-if(grepl('^copykat',method)) split(ix,ceiling(seq_along(ix)/20)) else if(method=='scevan') split(ix,ceiling(seq_along(ix)/100)) else as.list(ix)
    for(g in groups) {
      n<-n+1L;ans[[n]]<-if(length(g)>1) apply(v[g,,drop=FALSE],2,median) else v[g,,drop=TRUE]
      meta[[n]]<-data.frame(chr=cc,source_bins=length(g))
    }
  }
  z<-do.call(rbind,ans);colnames(z)<-ids
  payload<-data.table(do.call(rbind,meta),z)
  setnames(payload,c('chr','source_bins',colnames(z)))
  fwrite(payload,file.path(out,paste0(s,'_',method,'_profiles.tsv.gz')),sep='\t')
  cat('PLOT_DATA_READY',s,method,nrow(z),ncol(z),'\n')
}
cat('CNV_PLOT_DATA_EXPORTED\n')
