options(stringsAsFactors=FALSE)
root <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase51_matrix_lineage_attribution_20261002"
readtab <- function(n) read.delim(file.path(root,"tables",n),check.names=FALSE,na.strings=c("NA",""),quote="")
s <- read.delim(gzfile(file.path(root,"tables/sample_group_expression.tsv.gz")),check.names=FALSE,na.strings=NULL,quote="")
for(n in c("counts","library","cells","detected_cells"))s[[n]]<-as.numeric(s[[n]])
key <- function(x,cols) do.call(paste,c(x[cols],sep="\034"))
splitids <- split(seq_len(nrow(s)),key(s,c("dataset","patient_id","site","group","gene")))
site <- do.call(rbind,lapply(splitids,function(ii){
 t<-s[ii,];r<-t[1,c("dataset","patient_id","site","group","gene","group_scheme")]
 valid<-all(is.finite(t$counts[t$cells>0]))
 n<-sum(t$cells);lib<-sum(t$library)
 r$counts<-if(all(is.finite(t$counts)))sum(t$counts) else NA_real_
 r$library<-lib;r$cells<-n
 r$log_expression<-if(n>0 && lib>0 && valid)log2(sum(t$counts,na.rm=TRUE)/lib*1e4+.1) else NA_real_
 r
}))
pat <- readtab("patient_matched_effects.tsv")
calc <- rep(NA_real_,nrow(pat))
for(i in seq_len(nrow(pat))){
 r<-pat[i,];d<-site[site$dataset==r$dataset & site$patient_id==r$patient_id & site$gene==r$gene,]
 a<-d[d$group==r$target,c("site","log_expression")];b<-d[d$group==r$receiver,c("site","log_expression")]
 j<-merge(a,b,by="site",suffixes=c("_target","_receiver"))
 v<-j$log_expression_target-j$log_expression_receiver;v<-v[is.finite(v)]
 if(length(v))calc[i]<-mean(v)
}
stopifnot(identical(is.na(calc),is.na(pat$effect)),max(abs(calc-pat$effect),na.rm=TRUE)<1e-12)
pq<-data.frame(pat[,c("dataset","patient_id","target","receiver","gene")],R_effect=calc,Python_effect=pat$effect,abs_error=abs(calc-pat$effect))
write.table(pq,file.path(root,"tables/independent_R_patient_effect_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
st<-readtab("lineage_statistics.tsv");p2<-rep(NA_real_,nrow(st));m2<-rep(NA_real_,nrow(st))
for(i in seq_len(nrow(st))){
 r<-st[i,];ii<-pat$dataset==r$dataset & pat$target==r$target & pat$receiver==r$receiver & pat$gene==r$gene
 d<-calc[ii];d<-d[is.finite(d)];stopifnot(length(d)==r$N)
 if(length(d))m2[i]<-mean(d)
 if(length(d)>=2 && any(d!=0)){
  stopifnot(length(d)<=20)
  sums<-0
  for(v in d)sums<-c(sums+v,sums-v)
  p2[i]<-mean(abs(sums/length(d))>=abs(mean(d))-1e-12)
 }
}
q2<-rep(NA_real_,nrow(st))
for(f in unique(st$family)){
 ii<-which(st$family==f);stopifnot(length(ii)==27)
 q2[ii]<-p.adjust(p2[ii],method="BH",n=length(ii))
}
stopifnot(identical(is.na(p2),is.na(st$P)),identical(is.na(q2),is.na(st$FDR)))
qa<-data.frame(st[,c("dataset","target","receiver","gene")],mean_error=abs(m2-st$mean_effect),P_error=abs(p2-st$P),FDR_error=abs(q2-st$FDR))
stopifnot(max(qa$mean_error,na.rm=TRUE)<1e-12,max(qa$P_error,na.rm=TRUE)<1e-12,max(qa$FDR_error,na.rm=TRUE)<1e-12)
write.table(qa,file.path(root,"tables/independent_R_statistics_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
# Independently rebuild captured transcript shares from disjoint groups.
part<-site[site$group_scheme=="exhaustive_partition",]
totals<-aggregate(part$counts,part[,c("dataset","patient_id","site","gene")],sum)
names(totals)[5]<-"all_gene_counts"
sj<-merge(site,totals,by=c("dataset","patient_id","site","gene"))
sj$share<-ifelse(sj$all_gene_counts>0,sj$counts/sj$all_gene_counts,NA_real_)
ids<-split(seq_len(nrow(sj)),key(sj,c("dataset","patient_id","group","gene")))
pp<-do.call(rbind,lapply(ids,function(ii){
 t<-sj[ii,];r<-t[1,c("dataset","patient_id","group","gene")];d<-t$share[is.finite(t$share)]
 r$R_share<-if(length(d))mean(d) else NA_real_;r
}))
py<-readtab("patient_expression_and_shares.tsv");j<-merge(pp,py,by=c("dataset","patient_id","group","gene"))
stopifnot(nrow(j)==nrow(py),identical(is.na(j$R_share),is.na(j$gene_capture_share)))
j$abs_error<-abs(j$R_share-j$gene_capture_share)
stopifnot(max(j$abs_error,na.rm=TRUE)<1e-12)
write.table(j[,c("dataset","patient_id","group","gene","R_share","gene_capture_share","abs_error")],file.path(root,"tables/independent_R_RNA_share_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
out<-data.frame(check=c("patient_effect_rows","statistics_rows","RNA_share_rows","max_patient_effect_error","max_P_error","max_FDR_error","max_RNA_share_error"),value=c(nrow(pq),nrow(qa),nrow(j),max(pq$abs_error,na.rm=TRUE),max(qa$P_error,na.rm=TRUE),max(qa$FDR_error,na.rm=TRUE),max(j$abs_error,na.rm=TRUE)))
write.table(out,file.path(root,"delivery/INDEPENDENT_R_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines("SUCCESS",file.path(root,"INDEPENDENT_R_QA_COMPLETE"))
print(as.data.frame(out))
