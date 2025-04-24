# ================================================
# Seurat preprocessing script for HGSOC scRNA-seq
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
cellranger_dir <- file.path(base_dir, "output/cellranger")
metadata_path <- file.path(base_dir, "metadata/GSE184880_metadata.csv")
output_dir <- file.path(base_dir, "output/seurat")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load metadata
metadata <- read_csv(metadata_path)
print("Metadata loaded:")
print(head(metadata))

# Group runs by Sample Name
sample_groups <- split(metadata$Run, metadata$`Sample Name`)
print(paste("Found", length(sample_groups), "unique GSM samples"))

# Create list to store sample Seurat objects
seurat_samples <- list()

# Process each GSM sample by merging its technical replicates
for (sample_name in names(sample_groups)) {
  
  # Get SRR run IDs for this sample
  run_ids <- sample_groups[[sample_name]]
  print(paste("Processing", sample_name, "with", length(run_ids), "sequencing runs"))
  
  # List to store run-level Seurat objects
  run_objects <- list()
  
  # Process each run
  for (i in seq_along(run_ids)) {
    run_id <- run_ids[i]
    
    # Path to CellRanger output for this run
    matrix_dir <- file.path(cellranger_dir, run_id, "outs", "filtered_feature_bc_matrix")
    
    # Check if CellRanger output exists
    if (!dir.exists(matrix_dir)) {
      print(paste("Warning: CellRanger output for", run_id, "not found at", matrix_dir))
      next
    }
    
    # Load data
    print(paste("Loading", run_id))
    data <- Read10X(data.dir = matrix_dir)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = data,
      project = sample_name,
      min.cells = 3,
      min.features = 200
    )
    
    # Add metadata
    run_meta <- metadata[metadata$Run == run_id, ]
    seurat_obj$run_id <- run_id
    seurat_obj$sample_name <- sample_name
    seurat_obj$age <- run_meta$age[1]
    seurat_obj$pathology <- run_meta$pathology[1]
    seurat_obj$tissue_type <- run_meta$tissue_type[1]
    seurat_obj$tumor_stage <- run_meta$tumor_stage[1]
    
    # Calculate mitochondrial percentage
    seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    # Store in list with unique run identifier
    run_objects[[i]] <- seurat_obj
  }
  
  # If no data found for any run of this sample, skip
  if (length(run_objects) == 0) {
    print(paste("No data found for", sample_name, "- skipping"))
    next
  }
  
  # If only one run, use it directly, otherwise merge
  if (length(run_objects) == 1) {
    merged_sample <- run_objects[[1]]
    print(paste("Single run for", sample_name, "with", ncol(merged_sample), "cells"))
  } else {
    # Merge all runs for this sample
    print(paste("Merging", length(run_objects), "runs for", sample_name))
    
    # Add run ID to cell names to avoid conflicts when merging
    for (i in seq_along(run_objects)) {
      run_objects[[i]] <- RenameCells(run_objects[[i]], add.cell.id = run_ids[i])
    }
    
    # Merge runs
    merged_sample <- merge(run_objects[[1]], y = run_objects[2:length(run_objects)], 
                          add.cell.ids = NULL, project = sample_name)
    
    print(paste("Merged dataset for", sample_name, "has", ncol(merged_sample), "cells"))
  }
  
  # Store merged sample
  seurat_samples[[sample_name]] <- merged_sample
  
  # Plot QC metrics for this sample
  pdf(file.path(output_dir, paste0(sample_name, "_qc_plots.pdf")), width = 12, height = 8)
  
  p1 <- VlnPlot(merged_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 <- FeatureScatter(merged_sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p3 <- FeatureScatter(merged_sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  print(p1)
  print(p2)
  print(p3)
  
  dev.off()
  
  # Save this sample's Seurat object
  saveRDS(merged_sample, file.path(output_dir, paste0(sample_name, "_raw.rds")))
}

print(paste("Processed", length(seurat_samples), "samples"))
