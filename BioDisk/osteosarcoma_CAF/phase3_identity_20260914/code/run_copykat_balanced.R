# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' run_copykat_balanced.R
# RUN_COMMAND: Rscript --vanilla run_copykat_balanced.R <infercnv|copykat> <sample_short>
# EXPECTED_SUCCESS_MARKER: CNV_SAMPLE_COMPLETE
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L,args[1] %in% c('infercnv','copykat'))
method <- args[1]; sample_id <- args[2]
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914'
.libPaths(c(file.path(root,'Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
suppressPackageStartupMessages({library(Matrix);library(data.table);library(jsonlite)})
options(scipen=100,stringsAsFactors=FALSE)
set.seed(270231)
input <- file.path(root,'inputs',sample_id)
stopifnot(method=='copykat')
out <- file.path(root,'copykat_balanced',sample_id)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'COMPLETE'))) {cat('CNV_SAMPLE_CACHED',method,sample_id,'\n');quit(status=0)}
started <- format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')
status <- function(state,extra=list()) write_json(c(list(method=method,sample=sample_id,state=state,started=started,updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),extra),file.path(out,'status.json'),auto_unbox=TRUE,pretty=TRUE)
status('running')
tryCatch({
 stopifnot(file.exists(file.path(root,'resources','coordinate_audit.tsv')))
 cells <- as.data.frame(fread(file.path(input,'cells.tsv')))
 genes <- as.data.frame(fread(file.path(input,'genes.tsv')))
 mat <- as(readMM(gzfile(file.path(input,'counts.mtx.gz'))),'CsparseMatrix')
 dimnames(mat) <- list(genes$ensembl_id,cells$tool_id)
 stopifnot(ncol(mat)==nrow(cells),nrow(mat)==nrow(genes),!anyDuplicated(rownames(mat)),!anyDuplicated(colnames(mat)))
 stopifnot(all(is.finite(mat@x)),all(mat@x>=0),all(mat@x==floor(mat@x)),all(abs(colSums(mat)-cells$n_counts)<1e-6))
 refs_original <- cells$tool_id[grepl('^ref_',cells$cnv_role)]
 refs_immune <- cells$tool_id[cells$cnv_role=='ref_immune']
 refs_endo <- cells$tool_id[cells$cnv_role=='ref_endo']
 n_bal <- min(length(refs_immune),length(refs_endo))
 stopifnot(n_bal>=20)
 set.seed(270231)
 refs <- c(sample(refs_immune,n_bal),sample(refs_endo,n_bal))
 unused <- cells$tool_id %in% setdiff(refs_original,refs)
 cells$cnv_role[unused] <- 'check_immune'
 write.table(data.frame(tool_id=refs),file.path(out,'balanced_reference_cells.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 hold <- cells$tool_id[grepl('^holdout_',cells$cnv_role)]
 stopifnot(length(refs)>=30)
 result <- cells[,c('tool_id','cell_id','sample','short','patient_id','barcode','cluster','phase2_label','lineage_phase3','cnv_role','n_counts','n_genes')]
 result$method <- method;result$call <- 'not_evaluable';result$signal_strength <- NA_real_
 params <- list(method=method,sample=sample_id,seed=270231,reference='same_sample_immune_and_endothelial',reference_cells=length(refs),heldout_cells=length(hold),input_cells=ncol(mat),input_genes=nrow(mat),num_threads=4,coordinate_build='GRCh38',coordinate_annotation='Ensembl84',gene_id_type='Ensembl')
 if(method=='infercnv') {
   suppressPackageStartupMessages(library(infercnv))
   params$package_version <- as.character(packageVersion('infercnv'))
   params$cutoff <- .1;params$window_length <- 101;params$denoise <- TRUE;params$HMM <- FALSE
   # Full counts, gene-position based filtering by the tool; no HVG or integrated input.
   ic <- CreateInfercnvObject(raw_counts_matrix=mat,annotations_file=file.path(input,'annotations.tsv'),delim='\t',gene_order_file=file.path(root,'resources','infercnv_gene_order_GRCh38_84.tsv'),ref_group_names=unique(cells$cnv_role[grepl('^ref_',cells$cnv_role)]),chr_exclude=c('X','Y','MT','chrX','chrY','chrM'))
   ic <- infercnv::run(ic,cutoff=.1,min_cells_per_gene=3,out_dir=out,window_length=101,cluster_by_groups=TRUE,cluster_references=TRUE,denoise=TRUE,HMM=FALSE,analysis_mode='samples',num_threads=4,plot_steps=FALSE,inspect_subclusters=FALSE,resume_mode=TRUE,save_rds=FALSE,save_final_rds=TRUE,no_prelim_plot=TRUE,write_expr_matrix=FALSE,write_phylo=FALSE,output_format='png',png_res=150)
   stopifnot(nrow(ic@expr.data)>=2000)
   v <- log2(pmax(ic@expr.data,1e-3))
   score <- colMeans(abs(v));names(score)<-colnames(v)
   go <- ic@gene_order
   chrom <- as.character(go$chr)
   chrmed <- do.call(rbind,lapply(unique(chrom),function(cc) apply(v[chrom==cc,,drop=FALSE],2,median)))
   rownames(chrmed)<-unique(chrom);colnames(chrmed)<-colnames(v)
   # Independent heldout normal cells calibrate the descriptive signal cutoff.
   href <- intersect(hold,names(score));train <- intersect(refs,names(score))
   calibration <- if(length(href)>=20) href else train
   q99 <- unname(quantile(score[calibration],.99))
   robust <- median(score[calibration])+5*max(mad(score[calibration],constant=1),.005)
   threshold <- max(.05,q99,robust)
   broad <- colSums(abs(chrmed)>=.1)>=2 | apply(abs(chrmed),2,max)>=.25
   calls <- ifelse(score>threshold & broad,'elevated_signal',ifelse(score<=threshold,'low_signal','ambiguous_signal'))
   mapped <- match(result$tool_id,names(score));ok<-!is.na(mapped)
   result$signal_strength[ok]<-score[mapped[ok]];result$call[ok]<-calls[mapped[ok]]
   result$signal_cutoff <- threshold
   result$altered_chromosomes <- NA_integer_;result$altered_chromosomes[ok]<-colSums(abs(chrmed)>=.1)[mapped[ok]]
   fwrite(as.data.table(cbind(gene=rownames(v),ic@gene_order)),file.path(out,'gene_positions_used.tsv'),sep='\t')
   fwrite(data.table(chromosome=rownames(chrmed),chrmed),file.path(out,'chromosome_log2_medians.tsv'),sep='\t')
   # Coherent genomic windows retained as a lightweight profile for cross-method review.
   bins <- list();meta<-list();j<-0L
   for(cc in unique(chrom)) {
     ix<-which(chrom==cc);parts<-split(ix,ceiling(seq_along(ix)/100))
     for(ii in parts) {
       j<-j+1L;bins[[j]]<-apply(v[ii,,drop=FALSE],2,median)
       meta[[j]]<-data.frame(chr=cc,start=min(go$start[ii]),end=max(go$stop[ii]),genes=length(ii))
     }
   }
   bm<-do.call(rbind,bins);colnames(bm)<-colnames(v)
   saveRDS(list(bins=do.call(rbind,meta),matrix=bm),file.path(out,'window_profiles.rds'))
   params$mapped_genes <- nrow(v);params$signal_cutoff<-threshold
   params$calibration<-if(length(href)>=20)'heldout_normals' else 'training_normals_insufficient_holdout'
   params$signal_definition<-'mean_abs_log2_denoised_expression; elevation also requires >=2 chromosomes abs median>=0.1 or one>=0.25'
 } else {
   suppressPackageStartupMessages(library(copykat))
   params$package_version<-as.character(packageVersion('copykat'))
   params$reference_variant<-'equal_number_immune_and_endothelial; original heldout normals untouched'
   params$effective_internal_seed<-1234
   params$ngene_chr<-5;params$win_size<-25;params$KS_cut<-.1;params$LOW_DR<-.05;params$UP_DR<-.1
   params$genome_argument<-'hg20';params$coordinate_annotation<-'CopyKAT full.anno hg38, stable Ensembl-ID mapping'
   setwd(out)
   fulla<-get('full.anno',envir=as.environment('package:copykat'))
   amap<-fulla[fulla$ensembl_gene_id %in% rownames(mat),,drop=FALSE]
   amap<-amap[!duplicated(amap$hgnc_symbol),,drop=FALSE]
   ambiguous<-unique(amap$ensembl_gene_id[duplicated(amap$ensembl_gene_id)|duplicated(amap$ensembl_gene_id,fromLast=TRUE)])
   amap<-amap[!amap$ensembl_gene_id %in% ambiguous,,drop=FALSE]
   ckmat<-mat[as.character(amap$ensembl_gene_id),,drop=FALSE]
   stopifnot(nrow(ckmat)>5000)
   params$mapped_input_genes<-nrow(ckmat)
   params$annotation_ambiguous_ensembl_excluded<-length(ambiguous)
   fwrite(amap,file.path(out,'copykat_gene_mapping.tsv'),sep='\t')
   fwrite(data.table(ensembl_id=setdiff(rownames(mat),rownames(ckmat))),file.path(out,'copykat_unmapped_or_duplicate_genes.tsv'),sep='\t')
   ck <- copykat::copykat(rawmat=as.matrix(ckmat),id.type='E',ngene.chr=5,win.size=25,KS.cut=.1,sam.name=sample_id,distance='euclidean',norm.cell.names=refs,output.seg=FALSE,plot.genes=FALSE,genome='hg20',n.cores=4,LOW.DR=.05,UP.DR=.1,cell.line='no')
   saveRDS(ck,file.path(out,'copykat_result.rds'))
   pred<-as.data.frame(ck$prediction);stopifnot(all(c('cell.names','copykat.pred') %in% names(pred)))
   mm<-match(result$tool_id,pred$cell.names);ok<-!is.na(mm)
   result$call[ok]<-as.character(pred$copykat.pred[mm[ok]])
   result$call[!result$call %in% c('aneuploid','diploid')]<-'not_evaluable'
   cna<-as.data.frame(ck$CNAmat,check.names=FALSE)
   cc<-intersect(result$tool_id,names(cna));stopifnot(length(cc)>0)
   value<-as.matrix(cna[,cc,drop=FALSE]);stopifnot(all(is.finite(value)))
   result$signal_strength[match(cc,result$tool_id)]<-colMeans(abs(value))
   saveRDS(list(bins=cna[,setdiff(names(cna),cc),drop=FALSE],matrix=value),file.path(out,'window_profiles.rds'))
   params$bin_count<-nrow(value)
 }
 fwrite(result,file.path(out,'cell_results.tsv'),sep='\t',na='NA')
 write_json(params,file.path(out,'parameters.json'),auto_unbox=TRUE,pretty=TRUE)
 sink(file.path(out,'sessionInfo.txt'));print(sessionInfo());sink()
 status('complete',list(cells=nrow(result),calls=as.list(table(result$call))))
 writeLines(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),file.path(out,'COMPLETE'))
 cat('CNV_SAMPLE_COMPLETE',method,sample_id,'\n')
},error=function(e){
 msg<-conditionMessage(e);status('failed',list(error=msg));writeLines(msg,file.path(out,'FAILED'));cat('CNV_SAMPLE_FAILED',method,sample_id,msg,'\n');quit(status=1)
})
