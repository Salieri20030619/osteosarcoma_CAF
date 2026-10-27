library(Seurat)
library(copykat)
library(tidyverse)
library(dplyr)

# Load merged data
merge_obj <- readRDS("/mnt/d12b/NSCLC/bone/data/01.integrated/integrated/os_integrated_seuv5_rpca_v3_annotated.rds")

# Extract sample to list
# Assuming your merged Seurat object is called 'merged_seurat'
sample_list <- SplitObject(merge_obj, split.by = "seq_folder")

# To check the names of the samples extracted
names(sample_list)

# Load the seurat list containing all samples
# sample_list <- readRDS("/mnt/d12b/SC_BONE_RAW/data/sample_merged_processed/osteosarcoma_list_from_seurat_integration.rds")

# Define the output directory for saving the processed Seurat objects
# output_dir <- "/mnt/d12b/NSCLC/bone/data/05.CopyKat/copykat_v1"
# output_dir <- "/mnt/d12b/NSCLC/bone/data/05.CopyKat/copykat_v2"
# output_dir <- "/mnt/d12b/NSCLC/bone/data/05.CopyKat/copykat_v3"
output_dir <- "/mnt/d12b/NSCLC/bone/data/05.CopyKat/copykat_v5_manifest"

# # Increase mem
library(future)
#plan("multicore", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

# sample_list <- sample_list[!names(sample_list) %in% c("metastasis_2")] # No Metastasis_2
# names(sample_list)

# sample_list <- sample_list[c("metastasis_2")] # In case: We do separately for the metastasis_2
# names(sample_list)

# Loop for each sample in the list
for (sample_name in names(sample_list)) {
  cat("Processing sample:", sample_name, "\n")

  # Extract the Seurat object for the current sample
  seurat_obj <- sample_list[[sample_name]]

  # Join Layers 
  seurat_obj <- JoinLayers(seurat_obj)

  ## Prepare the Expression Matrix from the SCT assay counts
  # expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  expression_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "counts")


  # Define normal cells (e.g., endothelial cells and TILs) based on the cell_type metadata.
  # Adjust the cell types if needed.
  norm_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$rpca_cell_type %in% c("EC", "TILs"), ]) # Change to Myeloid and TILs when do for Metastasis_2
  #Select norm cells with rpca_cell_type == "EC" and nUMI > 2000
  norm_cells <- rownames(seurat_obj@meta.data[
    seurat_obj@meta.data$rpca_cell_type %in% c("EC", "TILs") & 
    seurat_obj@meta.data$nUMI > 2000 &
    seurat_obj@meta.data$mitoRatio < 0.02, 
  ])

  ## Run CopyKAT on the current sample
  copykat_results <- copykat(
    rawmat = expression_matrix,    # Expression matrix (genes x cells)
    sam.name = sample_name,
    distance = "euclidean",
    id.type = "S",                 # "S" for sample
    ngene.chr = 5,                 # Number of genes per chromosome
    win.size = 25,                 # Window size for smoothing
    KS.cut = 0.1,                  # Kolmogorov-Smirnov test cutoff
    output.seg = "FALSE",
    plot.genes = "TRUE",
    n.cores = 30,                  # Adjust number of cores as needed
    genome = "hg20",
    norm.cell.names = norm_cells
  )

  ## Extract the CopyKAT prediction results and prepare for merging with Seurat metadata
  pred_df <- data.frame(copykat_results$prediction)
  colnames(pred_df)[colnames(pred_df) == "cell.names"] <- "cells"

  # Add a column for cell barcodes if not already present
  seurat_obj@meta.data$cells <- rownames(seurat_obj@meta.data)

  # Merge CopyKAT results with the Seurat metadata by matching cell barcodes
  merged_meta <- merge(seurat_obj@meta.data, pred_df, by = "cells", all.x = TRUE)
  
  # Replace any "not.defined" calls with NA
  merged_meta$copykat.pred[merged_meta$copykat.pred == "not.defined"] <- NA
  rownames(merged_meta) <- merged_meta$cells

  # Update the Seurat object metadata with the merged information
  seurat_obj@meta.data <- merged_meta

  # ## Re-normalize the Seurat object using SCTransform
  # seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "mitoRatio", verbose = FALSE)
  # ## Dimensionality Reduction and Clustering
  # seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  # seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  # seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  # seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  
  ## Save the updated Seurat object with CNV annotations
  save_path <- file.path(output_dir, paste0(sample_name, "_annotated_filtered_copykat.rds"))
  saveRDS(seurat_obj, save_path)
  
  cat("Finished processing sample:", sample_name, "\n\n")
}