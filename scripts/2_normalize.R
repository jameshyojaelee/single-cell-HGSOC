# ================================================
# Seurat normalization and integration for HGSOC scRNA-seq
# ================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
metadata_path <- file.path(base_dir, "metadata/GSE184880_metadata.csv")
output_dir <- file.path(base_dir, "output/seurat")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load metadata to get sample info
metadata <- read_csv(metadata_path)
print("Metadata loaded")

# Find all raw RDS files in the output directory
rds_files <- list.files(path = output_dir, pattern = "_raw\\.rds$", full.names = TRUE)
print(paste("Found", length(rds_files), "raw RDS files to process:"))
print(rds_files)

# Initialize list to store sample Seurat objects
seurat_samples <- list()

# Load each raw RDS file
for (rds_file in rds_files) {
  # Extract sample name from filename
  sample_name <- basename(rds_file) %>% 
    str_replace("_raw\\.rds$", "")
  
  print(paste("Loading", sample_name, "from", rds_file))
  
  # Load the Seurat object
  seurat_obj <- readRDS(rds_file)
  
  # Store in the samples list
  seurat_samples[[sample_name]] <- seurat_obj
}

print(paste("Loaded", length(seurat_samples), "samples for processing"))

# Filtering and normalization for each sample
filtered_samples <- list()

for (sample_name in names(seurat_samples)) {
  print(paste("Filtering and normalizing", sample_name))
  
  # Get the sample
  sample <- seurat_samples[[sample_name]]
  
  # Filtering cells based on QC metrics
  sample <- subset(sample, subset = nFeature_RNA > 200 & percent.mt < 40)
  
  print(paste("After filtering:", ncol(sample), "cells remain"))
  
  # Normalize data
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  sample <- ScaleData(sample, features = rownames(sample))
  
  # Run PCA
  sample <- RunPCA(sample, features = VariableFeatures(object = sample))
  
  # Determine significant PCs
  pdf(file.path(output_dir, paste0(sample_name, "_elbow_plot.pdf")), width = 8, height = 6)
  print(ElbowPlot(sample, ndims = 50))
  dev.off()
  
  # adjust based on elbow plot
  n_pcs <- 20
  
  # UMAP and clustering
  sample <- FindNeighbors(sample, dims = 1:n_pcs)
  sample <- FindClusters(sample, resolution = 0.5)
  sample <- RunUMAP(sample, dims = 1:n_pcs)
  
  # Save plots
  pdf(file.path(output_dir, paste0(sample_name, "_umap_plots.pdf")), width = 12, height = 10)
  
  p1 <- DimPlot(sample, reduction = "umap", group.by = "seurat_clusters")
  print(p1)
  
  # If there are multiple runs in this sample, also show by run_id
  if (length(unique(sample$run_id)) > 1) {
    p2 <- DimPlot(sample, reduction = "umap", group.by = "run_id")
    print(p2)
  }
  
  # Save marker gene plots if you have known markers
  # Adjust these gene names for your specific cell types of interest
  epithelial_markers <- c("EPCAM", "KRT8", "KRT18")
  immune_markers <- c("PTPRC", "CD3E", "CD4", "CD8A", "CD14", "CD68")
  fibroblast_markers <- c("COL1A1", "DCN", "LUM")
  
  all_markers <- c(epithelial_markers, immune_markers, fibroblast_markers)
  all_markers <- all_markers[all_markers %in% rownames(sample)]
  
  if (length(all_markers) > 0) {
    print(FeaturePlot(sample, features = all_markers, ncol = 3))
    print(VlnPlot(sample, features = all_markers, ncol = 3))
  }
  
  dev.off()
  
  # Find markers for each cluster
  sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Save marker genes
  write.csv(sample.markers, file.path(output_dir, paste0(sample_name, "_markers.csv")))
  
  # Save processed object
  saveRDS(sample, file.path(output_dir, paste0(sample_name, "_processed.rds")))
  
  # Store for integration
  filtered_samples[[sample_name]] <- sample
}