suppressPackageStartupMessages({library(data.table)})
p<-'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917'
for (sample in c('GSE152048__BC22','GSE152048__BC10')) {
 x<-readRDS(file.path(p,'infercnv',sample,'window_profiles.rds'))
 fwrite(as.data.table(x$bins),file.path(p,'report_figures_20260917/data',paste0(sample,'_bins.tsv')),sep='\t')
 fwrite(data.table(tool_id=colnames(x$matrix),t(x$matrix)),file.path(p,'report_figures_20260917/data',paste0(sample,'_profiles.tsv.gz')),sep='\t',compress='gzip')
 cat('PROFILE_EXPORTED',sample,dim(x$matrix),'\n')
}
