# Osteosarcoma scRNA-seq Analysis
# OsteoCAR Atlas - Human Primary Tumors
# Data: Budhathoki et al. 2026 (https://doi.org/10.64898/2026.01.13.698472)

# Load libraries
library(qs)
library(Seurat)
library(ggplot2)

# Set working directory
setwd("~/Desktop/curtislab/introwork")

# Load Seurat object
sobj <- qread("patient_prim.qs")

# Explore metadata
colnames(sobj@meta.data)
table(sobj@meta.data$Ann_Level0)
table(sobj@meta.data$Ann_Level1)
table(sobj@meta.data$Ann_Level2)
table(sobj@meta.data$Ann_Level3)

# UMAP plots colored by annotation levels
p0 <- DimPlot(sobj, group.by = "Ann_Level0", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types - Level 0 (Tumor vs Host)")

p1 <- DimPlot(sobj, group.by = "Ann_Level1", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types - Level 1 (Broad)")

p2 <- DimPlot(sobj, group.by = "Ann_Level2", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types - Level 2 (Mid resolution)")

p3 <- DimPlot(sobj, group.by = "Ann_Level3", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types - Level 3 (Fine grained)")

# Save to PDF
pdf("patient_prim_umaps.pdf", width = 12, height = 10)
print(p0)
print(p1)
print(p2)
print(p3)
dev.off()