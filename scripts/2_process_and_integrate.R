# ================================================
# Seurat combined preprocessing, integration and analysis for HGSOC scRNA-seq
# ================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(future)

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

# STEP 1: Preprocessing each sample (minimal processing for integration)
filtered_samples <- list()

for (sample_name in names(seurat_samples)) {
  print(paste("Preprocessing", sample_name))
  
  # Get the sample
  sample <- seurat_samples[[sample_name]]
  
  # QC visualization before filtering
  pdf(file.path(output_dir, paste0(sample_name, "_qc_metrics.pdf")), width = 10, height = 8)
  p1 <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p3 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  # Filtering cells based on QC metrics
  sample <- subset(sample, subset = nFeature_RNA > 200 & percent.mt < 40)
  
  print(paste("After filtering:", ncol(sample), "cells remain"))
  
  # Normalize data
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  
  # Store preprocessed object for integration
  filtered_samples[[sample_name]] <- sample
  # save preprocessed object
  saveRDS(sample, file.path(output_dir, paste0(sample_name, "_preprocessed.rds")))
}

# STEP 2: Integration of all samples
if (length(filtered_samples) > 1) {
  print("Integrating all samples...")

  # Merge all samples into one object for IntegrateLayers
  merged <- merge(filtered_samples[[1]], y = filtered_samples[-1], add.cell.ids = names(filtered_samples))

  # Feature selection - benefits from parallelization
  plan("multicore", workers = 8)
  print("Finding variable features...")
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 1000)
  
  # Scale data - doesn't benefit as much from parallelization for most datasets
  plan("sequential")
  print("Scaling data...")
  merged <- ScaleData(merged)
  
  # PCA - benefits from parallelization
  plan("multicore", workers = 8)
  print("Running PCA...")
  merged <- RunPCA(merged, npcs = 30)

  # Integration - benefits from parallelization
  print("Integrating with CCAIntegration...")
  merged <- IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )

  # Join layers
  plan("sequential")
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
  print("Finding neighbors...")
  merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
  merged <- FindClusters(merged, resolution = 1)
  merged <- RunUMAP(merged, dims = 1:30, reduction = "integrated.cca")
  saveRDS(merged, file.path(output_dir, "integrated_seurat.rds"))

  print("Complete pipeline finished successfully!")
} else {
  print("ERROR: At least two samples are required for integration.")
}

# Ensure we reset to sequential at the end
plan("sequential")
print("Results saved to:")
print(output_dir)