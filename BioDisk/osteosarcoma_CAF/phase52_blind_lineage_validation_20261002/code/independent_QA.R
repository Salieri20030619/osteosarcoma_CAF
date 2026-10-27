options(stringsAsFactors=FALSE)
root <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase52_blind_lineage_validation_20261002"
readtab <- function(n) read.delim(file.path(root,"tables",n),check.names=FALSE,na.strings=c("NA",""),quote="")
cells <- read.delim(gzfile(file.path(root,"tables/cell_lineage_scores_and_predictions.tsv.gz")),check.names=FALSE,na.strings=c("NA",""),quote="")
raw <- read.delim(gzfile(file.path(root,"tables/raw_marker_counts_for_QA.tsv.gz")),check.names=FALSE,na.strings=c("NA",""),quote="")
panels <- jsonlite::fromJSON(file.path(root,"resources/LINEAGE_PANELS_FROZEN.json"))
genes <- c("COL1A1","LOX","MMP14");namespan <- names(panels)
stopifnot(length(intersect(genes,unlist(panels)))==0,identical(raw$cell_id,cells$cell_id),identical(raw$sample,cells$sample))
sets <- list(full=panels,A=lapply(panels,function(g)g[seq(1,length(g),by=2)]),B=lapply(panels,function(g)g[seq(2,length(g),by=2)]))
markercols <- sort(unique(unlist(panels)))
den <- raw$library_without_targets;den[den<=0] <- NA_real_
logmat <- log1p(as.matrix(raw[,markercols])/den*1e4)
getwinner <- function(v){
 safe<-v;safe[!is.finite(safe)] <- -Inf
 ix<-max.col(safe,ties.method="first");mx<-safe[cbind(seq_len(nrow(safe)),ix)]
 unique<-rowSums(abs(safe-mx)<=1e-12,na.rm=TRUE)==1 & is.finite(mx) & mx>0
 lab<-namespan[ix];lab[!unique]<-"Unresolved";lab
}
allscore <- list();labs <- list();scoresQA <- list()
for(split in names(sets)){
 v<-sapply(sets[[split]],function(g)rowMeans(logmat[,g,drop=FALSE],na.rm=TRUE))
 old<-as.matrix(cells[,paste0(split,"__",namespan)])
 stopifnot(identical(unname(is.na(v)),unname(is.na(old))))
 err<-max(abs(v-old),na.rm=TRUE);stopifnot(err<1e-12)
 lab<-getwinner(v);pycol<-switch(split,full="full_label",A="half_A_label",B="half_B_label")
 stopifnot(identical(lab,cells[[pycol]]))
 labs[[split]]<-lab
 scoresQA[[split]]<-data.frame(split=split,cells=nrow(cells),score_max_error=err,winner_mismatches=sum(lab!=cells[[pycol]]))
}
cons<-ifelse(labs$full==labs$A & labs$A==labs$B & labs$full!="Unresolved",labs$full,ifelse(labs$full=="Unresolved"|labs$A=="Unresolved"|labs$B=="Unresolved","Unresolved","Panel_discordant"))
stopifnot(identical(cons,cells$consensus_label))
scoreqa<-do.call(rbind,scoresQA);write.table(scoreqa,file.path(root,"tables/independent_R_score_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
tf<-function(v)tolower(as.character(v))%in%c("true","1")
records<-list();kk<-0
for(sm in unique(cells$sample)){
 ii<-which(cells$sample==sm);t<-cells[ii,]
 for(scheme in c("full_panel","split_panel_consensus")){
  lab<-if(scheme=="full_panel")labs$full[ii] else cons[ii]
  groupmasks<-list(query_Fibroblast=tf(t$query_all_unassigned)&lab=="Fibroblast",query_Osteogenic=tf(t$query_all_unassigned)&lab=="Osteogenic",supported_F3=tf(t$supported_F3))
  for(group in names(groupmasks)){
   mask<-groupmasks[[group]]
   for(g in genes){
    kk<-kk+1;records[[kk]]<-data.frame(scheme=scheme,sample=sm,patient_id=t$patient_id[1],site=t$site[1],group=group,gene=g,cells=sum(mask),counts=sum(t[[paste0(g,"_raw")]][mask]),library=sum(t$library[mask]))
   }
  }
 }
}
sg<-do.call(rbind,records)
site<-aggregate(sg[,c("cells","counts","library")],sg[,c("scheme","patient_id","site","group","gene")],sum)
site$log_expression<-ifelse(site$cells>0 & site$library>0,log2(site$counts/site$library*1e4+.1),NA_real_)
pat<-readtab("patient_target_gene_effects.tsv")
pairs<-list(query_Fibroblast_minus_query_Osteogenic=c("query_Fibroblast","query_Osteogenic"),query_Fibroblast_minus_supported_F3=c("query_Fibroblast","supported_F3"),query_Osteogenic_minus_supported_F3=c("query_Osteogenic","supported_F3"))
calc<-rep(NA_real_,nrow(pat))
for(i in seq_len(nrow(pat))){
 r<-pat[i,];d<-site[site$scheme==r$scheme&site$patient_id==r$patient_id&site$gene==r$gene,];pair<-pairs[[r$contrast]]
 a<-d[d$group==pair[1],c("site","log_expression")];b<-d[d$group==pair[2],c("site","log_expression")]
 j<-merge(a,b,by="site",suffixes=c("_L","_R"));v<-j$log_expression_L-j$log_expression_R;v<-v[is.finite(v)]
 if(length(v))calc[i]<-mean(v)
}
stopifnot(identical(is.na(calc),is.na(pat$effect)),max(abs(calc-pat$effect),na.rm=TRUE)<1e-12)
pq<-data.frame(pat[,c("scheme","patient_id","gene","contrast")],R_effect=calc,Python_effect=pat$effect,abs_error=abs(calc-pat$effect))
write.table(pq,file.path(root,"tables/independent_R_patient_effect_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
st<-readtab("target_gene_statistics.tsv");p2<-rep(NA_real_,nrow(st))
for(i in seq_len(nrow(st))){
 r<-st[i,];d<-calc[pat$scheme==r$scheme&pat$gene==r$gene&pat$contrast==r$contrast];d<-d[is.finite(d)]
 stopifnot(length(d)==r$N)
 if(length(d)>=2 & any(d!=0)){
  stopifnot(length(d)<=20);sums<-0
  for(v in d)sums<-c(sums+v,sums-v)
  p2[i]<-mean(abs(sums/length(d))>=abs(mean(d))-1e-12)
 }
}
q2<-rep(NA_real_,nrow(st))
for(f in unique(st$family)){ii<-which(st$family==f);stopifnot(length(ii)==9);q2[ii]<-p.adjust(p2[ii],"BH",n=length(ii))}
stopifnot(identical(is.na(p2),is.na(st$P)),identical(is.na(q2),is.na(st$FDR)),max(abs(p2-st$P),na.rm=TRUE)<1e-12,max(abs(q2-st$FDR),na.rm=TRUE)<1e-12)
stqa<-data.frame(st[,c("scheme","contrast","gene")],P_error=abs(p2-st$P),FDR_error=abs(q2-st$FDR))
write.table(stqa,file.path(root,"tables/independent_R_statistics_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
out<-data.frame(check=c("score_cells","panel_score_schemes","winner_mismatches","consensus_mismatches","max_score_error","patient_effect_rows","statistics_rows","max_patient_effect_error","max_P_error","max_FDR_error"),value=c(nrow(cells),3,0,0,max(scoreqa$score_max_error),nrow(pat),nrow(st),max(pq$abs_error,na.rm=TRUE),max(stqa$P_error,na.rm=TRUE),max(stqa$FDR_error,na.rm=TRUE)))
write.table(out,file.path(root,"delivery/INDEPENDENT_R_QA.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines("SUCCESS",file.path(root,"INDEPENDENT_R_QA_COMPLETE"));print(as.data.frame(out))
