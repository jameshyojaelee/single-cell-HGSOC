# ===============================================================
# Differential Expression Analysis using Seurat FindMarkers
# ===============================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(future)

# Increase the memory limit for future (important for large datasets)
# Adjust this value based on your system's memory
options(future.globals.maxSize = 200 * 1024^3) # 200GB limit, adjust as needed

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/seurat")
output_dir <- file.path(base_dir, "output/diff_expr")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the path to the integrated Seurat object
integrated_object_path <- file.path(input_dir, "integrated_seurat.rds")

# Check if the integrated object file exists
if (!file.exists(integrated_object_path)) {
  stop(paste("Integrated Seurat object not found at:", integrated_object_path))
}

# Load the integrated Seurat object
print(paste("Loading integrated Seurat object from:", integrated_object_path))
integrated_seurat <- readRDS(integrated_object_path)
print("Object loaded successfully.")

# Check for the 'tumor_stage' metadata column
if (!"tumor_stage" %in% colnames(integrated_seurat@meta.data)) {
    stop("Metadata column 'tumor_stage' not found in the Seurat object. Please check your metadata.")
}

# Print unique values in 'tumor_stage' for verification
print("Unique values in 'tumor_stage' column:")
print(unique(integrated_seurat@meta.data$tumor_stage))

# Ensure "Normal" exists in the tumor_stage column
if (!"Normal" %in% integrated_seurat@meta.data$tumor_stage) {
    stop("'Normal' group not found in the 'tumor_stage' column. Cannot perform comparison.")
}

# Set the cell identities to the 'tumor_stage' column for comparison
print("Setting cell identities based on 'tumor_stage'...")
Idents(integrated_seurat) <- "tumor_stage"

# Perform differential expression analysis: "Normal" vs. All Others
# Using Wilcoxon Rank Sum test
print("Running FindMarkers (Wilcoxon test) comparing 'Normal' vs. all other stages...")
plan("multicore", workers = 12) # Use multiple cores if available

# FindMarkers might take some time depending on the dataset size
markers <- FindMarkers(integrated_seurat, 
                       ident.1 = "Normal", 
                       test.use = "MAST", 
                       logfc.threshold = 0.25, # Default logfc threshold
                       min.pct = 0.1,          # Default min.pct threshold
                       verbose = TRUE)

plan("sequential") # Reset plan after parallel processing

# Check if markers were found
if (nrow(markers) == 0) {
  print("No significant markers found between 'Normal' and other stages with the current thresholds.")
} else {
  print(paste("Found", nrow(markers), "differential expression markers."))
  
  # Add gene symbols to the results (assuming they are in the row names)
  markers$gene <- rownames(markers)
  
  # Reorder columns for better readability (optional)
  markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  
  # Save the markers to a CSV file
  output_file_path <- file.path(output_dir, "diff_expr_Normal_vs_Others_wilcox.csv")
  print(paste("Saving differential expression results to:", output_file_path))
  write.csv(markers, file = output_file_path, row.names = FALSE)
  print("Results saved successfully.")
}

print("Differential expression analysis complete!")
print(paste("Output directory:", output_dir))
