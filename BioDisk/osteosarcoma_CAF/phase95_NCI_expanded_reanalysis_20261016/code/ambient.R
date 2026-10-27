args<-commandArgs(TRUE);d<-args[1]
suppressPackageStartupMessages({library(Matrix);library(data.table);library(SingleCellExperiment);library(celda);library(jsonlite)})
c<-as.data.frame(fread(file.path(d,'cells.tsv')));g<-as.data.frame(fread(file.path(d,'genes.tsv')))
x<-as(readMM(gzfile(file.path(d,'counts.mtx.gz'))),'CsparseMatrix');dimnames(x)<-list(g$feature_id,c$tool_id)
set.seed(20261016);sce<-SingleCellExperiment(list(counts=x));sce<-decontX(sce,z=as.integer(factor(c$cluster)),verbose=TRUE)
y<-decontXcounts(sce);lib<-colSums(y);z<-y%*%Diagonal(x=10000/pmax(lib,1));z@x<-log1p(z@x)
panels<-fromJSON('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase27_identity_review_20260926/resources/lineage_panels.json')
out<-data.frame(cell_id=c$cell_id,decontX_contamination=sce$decontX_contamination)
for(name in names(panels)){
 ix<-which(g$symbol%in%panels[[name]])
 out[[paste0('decontX_phase27_',name)]]<-if(length(ix)) as.numeric(colMeans(z[ix,,drop=FALSE])) else NA_real_
 out[[paste0('decontX_detected_',name)]]<-if(length(ix)) as.numeric(colSums(y[ix,,drop=FALSE]>=1)) else 0
}
scores<-as.matrix(out[,paste0('decontX_phase27_',names(panels))]);winner<-names(panels)[max.col(scores,ties.method='first')];mx<-apply(scores,1,max,na.rm=TRUE);winner[mx<=0|rowSums(scores==mx,na.rm=TRUE)!=1]<-'Unresolved';out$decontX_lineage_winner<-winner
fwrite(out,file.path(d,'ambient_evidence.tsv'),sep='\t');saveRDS(list(counts=y,contamination=sce$decontX_contamination,cluster=c$cluster),file.path(d,'decontX_result.rds'))
writeLines('AMBIENT_COMPLETE',file.path(d,'AMBIENT_COMPLETE'))
