.libPaths(c("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib", .libPaths()))

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
base <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF"
out <- file.path(base,"patient_expansion_download_investigation_20261013_v1")
.libPaths(c(file.path(out,"resources/Rlib"),.libPaths()))
stopifnot(requireNamespace("qs2",quietly=TRUE))
lineage <- fromJSON(file.path(base,"phase27_identity_review_20260926/resources/lineage_panels.json"))
programs <- fromJSON(file.path(base,"phase86_functional_bone_lung_20261012/resources/PROGRAMS.json"))
candidates <- fromJSON(file.path(base,"CAF_ALL14_RESEARCH_CURRENT.json"))$retained_genes
panels <- c(setNames(lineage,paste0("lineage__",names(lineage))),setNames(programs,paste0("program__",names(programs))))
checks <- list()
for (kind in c("patient_prim")) {
 cat("READ",kind,"\n");flush.console()
 obj <- qs2::qs_read(file.path(out,"downloads/OsteoCAR",paste0(kind,".qs2")),nthreads=2)
 stopifnot(inherits(obj,"Seurat"),!anyDuplicated(colnames(obj)))
 md <- obj@meta.data
 if("cell_id" %in% names(md)) md$author_cell_id_original <- md$cell_id
 md$cell_id <- rownames(md)
 for (n in names(md)) if(is.list(md[[n]])) md[[n]] <- vapply(md[[n]],function(z)paste(z,collapse=";"),character(1))
 con <- gzfile(file.path(out,"tables/OsteoCAR",paste0(kind,"_author_metadata.tsv.gz")),"wt")
 write.table(md,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(con)
 cat("METADATA",kind,nrow(md),paste(names(md),collapse="|"),"\n");flush.console()
 anns <- list()
 for (n in names(md)) if(is.character(md[[n]])||is.factor(md[[n]])) {
  v <- sort(table(md[[n]],useNA="ifany"),decreasing=TRUE)
  if(length(v)<=250) anns[[n]] <- as.list(v)
 }
 assay_names <- names(obj@assays)
 desc <- list(kind=kind,cells=ncol(obj),genes=nrow(obj),assays=assay_names,metadata_columns=names(md),category_counts=anns)
 raw_assay <- if("RNA" %in% assay_names) "RNA" else NA_character_
 if(!is.na(raw_assay)) {
  aa <- obj[[raw_assay]]
  layer_names <- if(inherits(aa,"Assay5")) SeuratObject::Layers(aa,search="^counts") else "counts"
  desc$raw_count_layers <- layer_names
  exports <- list();coverages <- list()
  for (layer in layer_names) {
   x <- if(inherits(aa,"Assay5")) SeuratObject::LayerData(aa,layer=layer) else SeuratObject::GetAssayData(obj,assay=raw_assay,layer="counts")
   stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(x@x>=0),all(is.finite(x@x)))
   lib <- Matrix::colSums(x);ng <- Matrix::colSums(x>0)
   y <- x %*% Matrix::Diagonal(x=10000/pmax(lib,1));y@x <- log1p(y@x)
   scores <- data.frame(cell_id=colnames(x),full_library_counts=as.numeric(lib),detected_features=as.numeric(ng),stringsAsFactors=FALSE)
   for (name in names(panels)) {
    genes <- panels[[name]];hit <- intersect(genes,rownames(x))
    scores[[name]] <- if(length(hit)) as.numeric(Matrix::colMeans(y[hit,,drop=FALSE])) else rep(NA_real_,ncol(x))
    coverages[[length(coverages)+1]] <- data.frame(layer=layer,program=name,defined=length(genes),measured=length(hit),fraction=length(hit)/length(genes),missing=paste(setdiff(genes,hit),collapse=";"))
   }
   for(g in candidates) scores[[paste0("gene_log1pCP10k__",g)]] <- if(g %in% rownames(x)) as.numeric(y[g,]) else NA_real_
   exports[[length(exports)+1]] <- scores
   desc[[paste0("layer_",layer)]] <- list(cells=ncol(x),genes=nrow(x),nonzero=length(x@x),gene_examples=head(rownames(x),12),integer_counts=all(abs(x@x-round(x@x))<1e-8))
   rm(x,y,scores);gc(FALSE)
  }
  if(length(exports)) {
   scores <- do.call(rbind,exports);stopifnot(!anyDuplicated(scores$cell_id),setequal(scores$cell_id,md$cell_id))
   con <- gzfile(file.path(out,"tables/OsteoCAR",paste0(kind,"_cell_scores.tsv.gz")),"wt");write.table(scores,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(con)
   write.table(do.call(rbind,coverages),file.path(out,"tables/OsteoCAR",paste0(kind,"_program_coverage.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
   desc$score_export_complete <- TRUE
  }
 }
 write_json(desc,file.path(out,"resources",paste0("OsteoCAR_",kind,"_INVENTORY.json")),pretty=TRUE,auto_unbox=TRUE,na="null")
 checks[[kind]] <- list(cells=ncol(obj),score_export=isTRUE(desc$score_export_complete),metadata_export=TRUE)
 cat("OBJECT_INVESTIGATION_COMPLETE",kind,"\n");flush.console()
 rm(obj,md,aa,exports,coverages);gc(FALSE)
}
write_json(list(status="completed",objects=checks,original_objects_unchanged=TRUE),file.path(out,"resources/OSTEOCAR_PRIM_EXTRACTION_COMPLETE.json"),pretty=TRUE,auto_unbox=TRUE)
