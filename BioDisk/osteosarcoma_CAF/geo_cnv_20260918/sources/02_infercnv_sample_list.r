library(infercnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Load the list of Seurat objects
# seurat_list <- readRDS("/mnt/d12b/SC_BONE_RAW/data/sample_merged_processed/osteosarcoma_list_from_merge.rds")

osteosarcoma_integration <- readRDS("/mnt/d12b/SC_BONE_RAW/data/sample_merged_processed/subclust/os_integrated_full_subclust.rds")
# Split the object
seurat_list <- SplitObject(osteosarcoma_integration, split.by = "seq_folder")


# Define constants for inferCNV
gene_order_file <- "/mnt/d12b/SC_BONE_RAW/data/table/input_infercnv/gene_ordering_hg38_genecode_v27.txt"
default_ref_groups <- c("EC", "TILs")
base_out_dir <- "/mnt/d12b/SC_BONE_RAW/data/out_infercnv_single_sample/update_v3/"

# Create the base output directory if it doesn't exist
if (!dir.exists(base_out_dir)) {
  dir.create(base_out_dir, recursive = TRUE)
}

# Initialize a list to store the inferCNV objects for each sample
infercnv_results_list <- list()

# Loop over each sample in the seurat list
for (sample_name in names(seurat_list)) {
  cat("Processing sample:", sample_name, "\n")
  
  # Load the Seurat object for the current sample
  seurat_obj <- seurat_list[[sample_name]]
  
  # Filter out cell types with fewer than 10 cells
  cell_counts <- table(seurat_obj$os_cell_type)
  valid_cells <- names(cell_counts)[cell_counts >= 10]
  seurat_obj <- subset(seurat_obj, os_cell_type %in% valid_cells)
  
  # Determine which reference groups are available in this sample.
  available_cell_types <- unique(seurat_obj$os_cell_type)
  ref_groups_used <- default_ref_groups[default_ref_groups %in% available_cell_types]
  
  if (length(ref_groups_used) == 0) {
    cat("Warning: None of the reference groups (", 
        paste(default_ref_groups, collapse = ", "),
        ") found in sample", sample_name, ". Skipping sample.\n\n")
    next  # Skip this sample if no valid reference groups are available
  } else {
    cat("Using reference group(s):", paste(ref_groups_used, collapse = ", "), "\n")
  }
  
  # Set the default assay to "RNA" for raw counts
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Extract the raw count matrix from the "RNA" assay
  raw_counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
  
  # Prepare the annotations data frame.
  # It must have rownames that match the cell barcodes in the raw_counts matrix.
  annotations <- seurat_obj@meta.data %>% 
    dplyr::select(os_cell_type)
  annotations <- annotations[colnames(raw_counts), , drop = FALSE]
  
  # Create the inferCNV object using the available reference groups
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = raw_counts,
                                         annotations_file = annotations,
                                         delim = "\t",
                                         gene_order_file = gene_order_file,
                                         ref_group_names = ref_groups_used)
  
  # Define an output directory for the current sample (this directory will be auto-created by infercnv::run)
  sample_out_dir <- file.path(base_out_dir, sample_name)
  
  # Run inferCNV for the current sample
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,            # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = sample_out_dir,  # output directory for this sample
                                cluster_by_groups = TRUE,  # cluster cells by their annotation groups
                                denoise = TRUE,
                                resume_mode = TRUE,
                                num_threads = 86,          # adjust as needed
                                HMM = TRUE)
  
  # Store the resulting inferCNV object in the list
  infercnv_results_list[[sample_name]] <- infercnv_obj
  
  cat("Finished processing sample:", sample_name, "\n\n")
}

# Optionally, save the list of inferCNV objects for downstream analyses
saveRDS(infercnv_results_list, file = file.path(base_out_dir, "all_infercnv_results.rds"))
cat("All samples have been processed and inferCNV results saved.\n")