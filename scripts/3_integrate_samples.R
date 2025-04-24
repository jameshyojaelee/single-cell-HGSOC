# ================================================
# Seurat integration for HGSOC scRNA-seq - Step 3
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
output_dir <- file.path(base_dir, "output/seurat")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Find all processed RDS files in the output directory
rds_files <- list.files(path = output_dir, pattern = "_processed\\.rds$", full.names = TRUE)
print(paste("Found", length(rds_files), "processed RDS files to integrate:"))
print(rds_files)

# Initialize list to store sample Seurat objects
filtered_samples <- list()

# Load each processed RDS file
for (rds_file in rds_files) {
  # Extract sample name from filename
  sample_name <- basename(rds_file) %>% 
    str_replace("_processed\\.rds$", "")
  
  print(paste("Loading", sample_name, "from", rds_file))
  
  # Load the Seurat object
  seurat_obj <- readRDS(rds_file)
  
  # Store in the samples list
  filtered_samples[[sample_name]] <- seurat_obj
}

print(paste("Loaded", length(filtered_samples), "samples for integration"))

# Integration of all samples if there are multiple samples
if (length(filtered_samples) > 1) {
  print("Integrating all samples...")
  
  # Prepare for integration
  sample.list <- filtered_samples
  features <- SelectIntegrationFeatures(object.list = sample.list)
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = sample.list, 
                                   anchor.features = features)
  
  # Create integrated dataset
  integrated <- IntegrateData(anchorset = anchors)
  
  # Switch to integrated assay for downstream analysis
  DefaultAssay(integrated) <- "integrated"
  
  # Scale the integrated data
  integrated <- ScaleData(integrated, verbose = FALSE)
  
  # Run PCA
  integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
  
  # UMAP and clustering
  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
  integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
  integrated <- FindClusters(integrated, resolution = 0.5)
  
  # Save plots
  pdf(file.path(output_dir, "integrated_analysis.pdf"), width = 12, height = 10)
  
  p1 <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters")
  p2 <- DimPlot(integrated, reduction = "umap", group.by = "sample_name")
  p3 <- DimPlot(integrated, reduction = "umap", group.by = "tissue_type")
  p4 <- DimPlot(integrated, reduction = "umap", group.by = "tumor_stage")
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  # Check for batch effects
  p5 <- DimPlot(integrated, reduction = "umap", split.by = "sample_name")
  print(p5)
  
  dev.off()
  
  # Find markers for integrated clusters
  DefaultAssay(integrated) <- "RNA"
  integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(integrated.markers, file.path(output_dir, "integrated_markers.csv"))
  
  # Save integrated object
  saveRDS(integrated, file.path(output_dir, "integrated_seurat.rds"))
  
  print("Integration complete!")
} else {
  print("ERROR: At least two samples are required for integration.")
}

print("Integration complete! Results saved to:")
print(output_dir)