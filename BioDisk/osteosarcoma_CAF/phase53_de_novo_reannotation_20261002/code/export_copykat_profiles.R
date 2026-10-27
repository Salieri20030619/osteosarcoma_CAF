options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(rhdf5);library(data.table)})
root<-"/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002"
dir.create(file.path(root,"objects","cnv"),recursive=TRUE,showWarnings=FALSE)
manifest<-read.delim(file.path(root,"resources","COPYKAT_PROFILE_MANIFEST.tsv"),check.names=FALSE,na.strings="")
qa<-list();metrics<-list()
for(i in seq_len(nrow(manifest))){
 r<-manifest[i,];cells<-read.delim(r$cells_table)
 if(is.na(r$copykat_rds)||!nzchar(r$copykat_rds)){
  qa[[i]]<-data.frame(sample=r$sample,input_cells=nrow(cells),profile_cells=0,bins=0,state="profile_unavailable")
  next
 }
 message("COPYKAT_PROFILE ",r$sample)
 ck<-readRDS(r$copykat_rds);tab<-ck$CNAmat
 stopifnot(all(c("chrom","chrompos","abspos")%in%names(tab)))
 cnv<-as.matrix(tab[,setdiff(names(tab),c("chrom","chrompos","abspos")),drop=FALSE])
 j<-match(colnames(cnv),cells$tool_id);stopifnot(!anyNA(j),!anyDuplicated(j))
 chr<-as.character(tab$chrom);chr[chr=="23"]<-"X";chr[chr=="24"]<-"Y"
 start<-floor(as.numeric(tab$chrompos)/5e6)*5e6
 bins<-paste0("chr",chr,":",format(start,scientific=FALSE,trim=TRUE))
 b<-rowsum(cnv,group=bins,reorder=FALSE);b<-b/as.numeric(table(factor(bins,levels=rownames(b))))
 stopifnot(all(is.finite(b)))
 meta<-data.frame(bin_id=rownames(b),chr=sub(":.*","",rownames(b)),start=as.numeric(sub(".*:","",rownames(b))),end=as.numeric(sub(".*:","",rownames(b)))+5e6)
 f<-file.path(root,"objects","cnv",paste0(r$sample,"_copykat_profiles.h5"))
 if(file.exists(f))stop("Output exists: retain completed export")
 h5createFile(f);h5write(b,f,"profiles");h5write(cells$cell_id[j],f,"cell_id");h5write(meta$bin_id,f,"bin_id")
 write.table(meta,file.path(root,"tables",paste0(r$sample,"_copykat_bins.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
 metrics[[i]]<-data.frame(sample=r$sample,cell_id=cells$cell_id[j],mean_abs_native_signal=colMeans(abs(b)),mean_native_signal=colMeans(b),max_abs_binned_signal=apply(abs(b),2,max),binned_regions=nrow(b))
 qa[[i]]<-data.frame(sample=r$sample,input_cells=nrow(cells),profile_cells=ncol(b),bins=nrow(b),state="exported")
 rm(ck,tab,cnv,b);gc(FALSE)
}
q<-do.call(rbind,qa);write.table(q,file.path(root,"tables","copykat_profile_coverage.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
d<-do.call(rbind,metrics);write.table(d,gzfile(file.path(root,"tables","copykat_cell_continuous_metrics.tsv.gz")),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(nrow(q)==107,sum(q$state=="exported")==84)
writeLines("SUCCESS",file.path(root,"COPYKAT_PROFILES_COMPLETE"))
cat("COPYKAT_PROFILES_COMPLETE",nrow(d),"\n")
