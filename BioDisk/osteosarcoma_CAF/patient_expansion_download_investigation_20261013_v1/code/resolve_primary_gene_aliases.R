
.libPaths(c("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1/resources/Rlib",.libPaths()))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(jsonlite))
O <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1"
x <- qs2::qs_read(file.path(O,"downloads/OsteoCAR/patient_prim.qs2"),nthreads=2)
c <- SeuratObject::GetAssayData(x,assay="RNA",layer="counts")
scores_file <- file.path(O,"tables/OsteoCAR/patient_prim_cell_scores.tsv.gz")
s <- read.delim(scores_file,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(setequal(s$cell_id,colnames(c)))
ix <- match(s$cell_id,colnames(c));aliases <- c("PEF1-AS1"="AC114488.2","RGS2-AS1"="AL390957.1")
file.copy(scores_file,file.path(O,"resources/patient_prim_cell_scores_before_alias_resolution.tsv.gz"),overwrite=FALSE)
for(g in names(aliases)) {stopifnot(aliases[[g]] %in% rownames(c));s[[paste0("gene_log1pCP10k__",g)]] <- log1p(as.numeric(c[aliases[[g]],ix])/s$full_library_counts*10000)}
s$available_matrix_library_counts <- s$full_library_counts
con <- gzfile(scores_file,"wt");write.table(s,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(con)
candidates <- fromJSON(file.path(dirname(O),"CAF_ALL14_RESEARCH_CURRENT.json"))$retained_genes
map <- do.call(rbind,lapply(candidates,function(g){feature <- if(g %in% rownames(c)) g else if(g %in% names(aliases) && aliases[[g]] %in% rownames(c)) aliases[[g]] else NA_character_;data.frame(canonical_gene=g,object_feature=feature,assayed=!is.na(feature))}))
write.table(map,file.path(O,"tables/OsteoCAR/patient_prim_candidate_feature_map.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write_json(list(status="complete",gene_mapping=map,original_QS_unchanged=TRUE),file.path(O,"resources/PRIMARY_ALIAS_RESOLUTION_COMPLETE.json"),pretty=TRUE,auto_unbox=TRUE)
print(as.data.frame(map))
