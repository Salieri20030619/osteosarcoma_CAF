args<-commandArgs(TRUE);root<-args[1];out<-file.path(root,'result_review_20260919','tables')
s<-read.delim(file.path(root,'review','input_summary.tsv'));s<-s[s$cnv_eligible=='True',]
for(sm in s$sample){
 p<-readRDS(file.path(root,'copykat',sm,'window_profiles.rds'))
 if(sm==s$sample[1]){print(as.data.frame(head(p$bins)));print(dim(p$matrix))}
 stopifnot('chrom' %in% names(p$bins),nrow(p$bins)==nrow(p$matrix),!anyDuplicated(colnames(p$matrix)))
 chr<-sub('^chr','',as.character(p$bins$chrom));keep<-as.character(1:22)
 med<-t(sapply(keep,function(cc)apply(p$matrix[chr==cc,,drop=FALSE],2,median)))
 stopifnot(all(is.finite(med)))
 write.table(data.frame(chromosome=keep,med,check.names=FALSE),gzfile(file.path(out,paste0(sm,'_copy_chromosomes.tsv.gz'))),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(chromosome=keep,bins=sapply(keep,function(cc)sum(chr==cc))),file.path(out,paste0(sm,'_copy_coverage.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 cat('PROFILE_EXPORT_COMPLETE',sm,ncol(med),'\n');rm(p,med);gc()
}
cat('ALL_PROFILES_EXPORTED\n')
