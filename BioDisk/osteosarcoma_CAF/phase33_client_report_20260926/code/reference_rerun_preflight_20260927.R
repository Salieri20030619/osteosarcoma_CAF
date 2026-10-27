suppressPackageStartupMessages(library(infercnv))
O<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase33_client_report_20260926'
for(sm in c('SCPCS000614','SCPCS000610'))for(v in c('source_reference','lymphoid_reference')){
 x<-readRDS(file.path(O,'cnv_reference',sm,v,'run.final.infercnv_obj'))
 a<-read.delim(file.path(O,'inputs',sm,paste0(v,'_annotations.tsv')),header=FALSE)
 actual<-colnames(x@expr.data)[unlist(x@reference_grouped_cell_indices,use.names=FALSE)]
 expected<-a[a[[2]]=='reference',1]
 stopifnot(setequal(actual,expected),!anyDuplicated(actual),length(intersect(unlist(x@reference_grouped_cell_indices),unlist(x@observation_grouped_cell_indices)))==0)
 cat(sm,v,'REFERENCE_MEMBERS_EXACT',length(actual),'OBSERVED',sum(lengths(x@observation_grouped_cell_indices)),'\n')
}
