# ================================================
# Seurat integration of preprocessed HGSOC samples
# ================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(future)

# Increase the memory limit for future (important for large datasets)
options(future.globals.maxSize = 200 * 1024^3)

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
output_dir <- file.path(base_dir, "output/seurat")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Find all preprocessed RDS files
preprocessed_files <- list.files(path = output_dir, pattern = "_preprocessed\\.rds$", full.names = TRUE)
print(paste("Found", length(preprocessed_files), "preprocessed files to integrate:"))
print(preprocessed_files)

if (length(preprocessed_files) < 2) {
  stop("At least two preprocessed samples are required for integration.")
}

# Load preprocessed objects
preprocessed_samples <- list()
for (rds_file in preprocessed_files) {
  sample_name <- basename(rds_file) %>% 
    str_replace("_preprocessed\\.rds$", "")
  
  print(paste("Loading", sample_name))
  preprocessed_samples[[sample_name]] <- readRDS(rds_file)
}

print(paste("Loaded", length(preprocessed_samples), "samples for integration"))

# INTEGRATION WORKFLOW
print("Starting integration workflow...")


# Merge all samples
print("Merging samples...")
merged <- merge(preprocessed_samples[[1]], 
                y = preprocessed_samples[-1], 
                add.cell.ids = names(preprocessed_samples))

# Free memory
rm(preprocessed_samples)
gc()

# Feature selection
plan("multicore", workers = 8)
print("Finding variable features...")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 1000)

# Scale data
plan("sequential")
print("Scaling data...")
merged <- ScaleData(merged)

# Run PCA
plan("multicore", workers = 8)
print("Running PCA...")
merged <- RunPCA(merged, npcs = 30)

# Save pre-integration object as checkpoint
print("Saving pre-integration checkpoint...")
saveRDS(merged, file.path(output_dir, "pre_integration_checkpoint.rds"))

# Switch to sequential for integration
plan("sequential")
print("Integrating with CCAIntegration (this may take a while)...")
tryCatch({
  merged <- IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    orig.reduction = "pca", 
    new.reduction = "integrated.cca",
    verbose = TRUE
  )
  
  # Save post-integration checkpoint
  print("Saving post-integration checkpoint...")
  saveRDS(merged, file.path(output_dir, "post_integration_checkpoint.rds"))
  
}, error = function(e) {
  print(paste("Integration error:", e$message))
  print("Saving current state before error...")
  saveRDS(merged, file.path(output_dir, "integration_error_checkpoint.rds"))
  stop("Integration failed. See error message above.")
})

# Re-join layers
print("Joining layers...")
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

# Downstream analysis
print("Running downstream analysis...")
merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
merged <- FindClusters(merged, resolution = 1)
merged <- RunUMAP(merged, dims = 1:30, reduction = "integrated.cca")

# Save final object
print("Saving integrated object...")
saveRDS(merged, file.path(output_dir, "integrated_seurat.rds"))

# Visualization
print("Creating visualizations...")
pdf(file.path(output_dir, "integrated_umap_plots.pdf"), width = 12, height = 10)

# Basic cluster plot
p1 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters")
print(p1)

# Sample visualization (if sample metadata exists)
if ("sample" %in% colnames(merged@meta.data) || "orig.ident" %in% colnames(merged@meta.data)) {
  sample_col <- ifelse("sample" %in% colnames(merged@meta.data), "sample", "orig.ident")
  p2 <- DimPlot(merged, reduction = "umap", group.by = sample_col)
  print(p2)
  
  # Split by sample
  p3 <- DimPlot(merged, reduction = "umap", split.by = sample_col)
  print(p3)
}

dev.off()

print("Integration and analysis complete!")
print("Results saved to:")
print(output_dir) 