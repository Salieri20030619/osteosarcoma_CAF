suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
p <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase97_normal_background_20261020"
x <- readRDS(file.path(p,"inputs/bone_msc_verified.rds"))
m <- x[[]]
a <- fromJSON(file.path(p,"resources/GENE_ALIAS_MAP.json"))
genes <- names(a$aliases)
stopifnot("cluster_anno_l2" %in% names(m))
labelcol <- "cluster_anno_l2"
masks <- list(as.character(m[[labelcol]])=="Fibro-MSC")
stopifnot(sum(masks[[1]])==1443)
donorcols <- names(m)[vapply(m,function(v)all(grepl("^H[0-9]+(_MACS)?$",as.character(v))),logical(1))]
if(!length(donorcols)) stop("Donor column requires explicit review")
donors <- lapply(donorcols,function(k)sub("_MACS$","",as.character(m[[k]])))
if(!all(vapply(donors,function(z)identical(z,donors[[1]]),logical(1)))) stop("Ambiguous donor columns")
donor <- donors[[1]]
labels <- as.character(m[[labelcol]])
raw <- GetAssayData(x,assay="RNA",layer="counts")
stopifnot(identical(colnames(raw),rownames(m)),all(raw@x>=0),all(raw@x==floor(raw@x)),all(is.finite(raw@x)))
lib <- Matrix::colSums(raw)
stopifnot("nCount_RNA" %in% names(m),all(lib==m$nCount_RNA))
coverage <- do.call(rbind,lapply(genes,function(g){
 matches <- intersect(rownames(raw),unique(c(g,unlist(a$aliases[[g]]),a$stable_ids[[g]])))
 if(length(matches)>1) stop(paste("Ambiguous gene",g,paste(matches,collapse=";")))
 data.frame(gene=g,reference_row=if(length(matches))matches else NA_character_,measured=length(matches)==1)
}))
write.table(coverage,file.path(p,"tables/bone_gene_coverage.tsv"),sep="\t",fileEncoding="UTF-8",quote=FALSE,row.names=FALSE,na="NA")
specs <- list("Fibro-MSC"=masks[[1]],"all_MSC"=grepl("MSC",labels,ignore.case=TRUE)&!grepl("RNA.*lo",labels,ignore.case=TRUE))
rows <- list();cellrows <- list()
for(spec in names(specs)){
 selected <- which(specs[[spec]])
 stopifnot(length(selected)>0)
 cellrows[[spec]] <- data.frame(specification=spec,cell_id=rownames(m)[selected],donor=donor[selected],published_label=labels[selected],library_UMIs=lib[selected])
 for(dd in sort(unique(donor[selected]))){
  ix <- selected[donor[selected]==dd]
  denom <- sum(lib[ix])
  stopifnot(denom>0)
  for(j in seq_len(nrow(coverage))){
   g <- coverage$gene[j];measured <- coverage$measured[j]
   count <- if(measured)sum(raw[coverage$reference_row[j],ix,drop=FALSE]) else NA_real_
   rows[[length(rows)+1]] <- data.frame(reference="bone",specification=spec,donor=dd,gene=g,raw_count=count,library_UMIs=denom,cells=length(ix),measured=measured,expression=if(measured)log2(count/denom*1e4+0.1) else NA_real_)
  }
 }
}
out <- do.call(rbind,rows)
write.table(out,file.path(p,"tables/bone_reference_pseudobulk.tsv"),sep="\t",fileEncoding="UTF-8",quote=FALSE,row.names=FALSE,na="NA")
write.table(do.call(rbind,cellrows),file.path(p,"tables/bone_selected_cell_metadata.tsv"),sep="\t",fileEncoding="UTF-8",quote=FALSE,row.names=FALSE)
write_json(list(annotation_column=labelcol,donor_column=donorcols[1],reference_gene_rows=nrow(raw),reference_cells=ncol(raw),raw_integer_nonnegative=TRUE,nCount_RNA_roundtrip=TRUE,published_labels=as.list(table(labels)),primary_donors=length(unique(out$donor[out$specification=="Fibro-MSC"])),primary_cells=sum(specs[["Fibro-MSC"]]),broad_cells=sum(specs[["all_MSC"]]),missing_genes=coverage$gene[!coverage$measured]),file.path(p,"resources/BONE_EXTRACTION_QA.json"),auto_unbox=TRUE,pretty=TRUE)
cat("BONE_EXTRACTION_DONE\n")
print(as.data.frame(table(out$specification)))
