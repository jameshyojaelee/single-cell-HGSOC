# =====================================================
# InferCNV Analysis for a Single Stage
# =====================================================
# This script runs InferCNV for one specific tumor stage,
# designed to be launched in parallel for multiple stages.

# --- Load required libraries ---
library(infercnv, lib="/gpfs/commons/home/jameslee/R/x86_64-pc-linux-gnu-library/4.4/")
library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_seurat_dir <- file.path(base_dir, "output/seurat")
output_infercnv_dir <- file.path(base_dir, "output/infercnv")

# Create output directory if it doesn't exist (InferCNV run will create stage-specific subdir)
if (!dir.exists(output_infercnv_dir)) dir.create(output_infercnv_dir, recursive = TRUE)

# Input Seurat object
annotated_object_path <- file.path(input_seurat_dir, "integrated_seurat_annotated.rds")

# Path to gene order file (3 columns: gene, chr, start, end, NO header)
gene_order_file_path <- "/gpfs/commons/home/jameslee/HGSOC/output/infercnv/hg38_gencode_v27.txt"

# Metadata column for cell type labels
cell_identity_col <- "sctype_classification"
# Metadata column for conditions/stages
condition_col <- "tumor_stage"
# Define reference stage (cells from this stage are used as normal references)
reference_stage <- "Normal"
# Define non-immune cell types based on your annotations
non_immune_cell_types <- c("Non-immune cells", "Endothelial") # Updated based on 4_celltype_annotation.R output

# InferCNV Parameters
infercnv_cutoff <- 0.1 # Default is 0.1; Use 1 for Smart-seq2, 0.1 for 10x Genomics
infercnv_window <- 101 # Default; sliding window size for smoothing
infercnv_threshold <- 0.1 # Default; minimum threshold for CNV calls
infercnv_denoise <- TRUE # Default; apply denoising
infercnv_HMM <- FALSE # Default; Set TRUE to run Hidden Markov Model prediction (slower, requires HMMiR)
num_threads <- 16 # Number of cores to use (per stage run)

# --- Get Stage Argument ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 11_infercnv_stage_runner.R <stage_name>", call. = FALSE)
}
stage <- args[1]
cat(paste("Received stage argument:", stage, "
"))


# --- Load Seurat Object and Check Files ---

if (!file.exists(annotated_object_path)) {
  stop(paste("Annotated Seurat object not found at:", annotated_object_path))
}
if (!file.exists(gene_order_file_path)) {
  stop(paste("Gene order file not found at:", gene_order_file_path,
             "
Please provide a valid path to a gene order file."))
}

cat(paste("Loading annotated Seurat object from:", annotated_object_path, "
"))
seurat_obj <- readRDS(annotated_object_path)
cat("Annotated object loaded successfully.
")

# Check required metadata
if (!cell_identity_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required cell identity column not found:", cell_identity_col))
}
if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required condition column not found:", condition_col))
}
if (!reference_stage %in% seurat_obj@meta.data[[condition_col]]) {
    stop(paste("Reference stage", reference_stage, "not found in", condition_col, "metadata."))
}
# Check if the provided stage exists in the data
all_stages <- unique(seurat_obj@meta.data[[condition_col]])
if (!stage %in% all_stages) {
    stop(paste("Provided stage '", stage, "' not found in the '", condition_col, "' metadata column. Available stages: ", paste(all_stages, collapse=", "), sep=""))
}
if (stage == reference_stage) {
    stop(paste("Provided stage '", stage, "' is the reference stage. This script processes tumor stages only.", sep=""))
}


# --- Run InferCNV for the Specified Stage ---

cat(paste("
======== Processing Stage:", stage, "========
"))

stage_output_dir <- file.path(output_infercnv_dir, stage)
if (!dir.exists(stage_output_dir)) dir.create(stage_output_dir)

# 1. Subset Seurat object: Current tumor stage + Reference stage
cat("  1. Subsetting Seurat object...
")
cells_to_keep <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[condition_col]] %in% c(stage, reference_stage), ])
seurat_subset <- subset(seurat_obj, cells = cells_to_keep)

# 2. Filter for non-immune cell types
cat("  2. Filtering for non-immune cell types...
")
non_immune_cells <- rownames(seurat_subset@meta.data[seurat_subset@meta.data[[cell_identity_col]] %in% non_immune_cell_types, ])
seurat_subset_filtered <- subset(seurat_subset, cells = non_immune_cells)

n_cells <- ncol(seurat_subset_filtered)
if (n_cells < 50) { # Arbitrary minimum cell count
    cat(paste("  Skipping stage", stage, ": Too few non-immune cells after filtering (", n_cells, ").
"))
    quit(status = 0) # Exit gracefully
}
cat(paste("     Kept", n_cells, "cells from types:", paste(unique(seurat_subset_filtered@meta.data[[cell_identity_col]]), collapse=", "), "
"))

# 3. Prepare Raw Counts Matrix
cat("  3. Preparing counts matrix...
")
counts_matrix <- GetAssayData(seurat_subset_filtered, assay = "RNA", slot = "counts")
# Ensure matrix is non-sparse (infercnv prefers dense)
if (is(counts_matrix, "sparseMatrix")) {
    counts_matrix <- as.matrix(counts_matrix)
}
# Remove genes with no expression in this subset
counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]
cat(paste("     Matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells
"))

# 4. Prepare Cell Annotations File
cat("  4. Preparing cell annotations file...
")
cell_annotations <- seurat_subset_filtered@meta.data %>%
  select(all_of(c(condition_col, cell_identity_col))) %>%
  # Create unique group names like "Normal_Epithelial" or "IC2_Fibroblast"
  mutate(group = paste(.data[[condition_col]], .data[[cell_identity_col]], sep="_")) %>%
  select(group) %>%
  rownames_to_column("cell_id")

annotations_file <- file.path(stage_output_dir, paste0("cell_annotations_", stage, ".txt"))
write.table(cell_annotations, file = annotations_file, sep = "	", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(paste("     Annotations file saved to:", annotations_file, "
"))

# 5. Define Reference Group Names
cat("  5. Defining reference groups...
")
reference_groups <- unique(cell_annotations$group[startsWith(cell_annotations$group, paste0(reference_stage, "_"))])
if (length(reference_groups) == 0) {
    cat(paste("  Skipping stage", stage, ": No reference cells found after filtering for non-immune types.
"))
    quit(status = 0) # Exit gracefully
}
cat(paste("     Using reference groups:", paste(reference_groups, collapse=", "), "
"))

# 6. Create InferCNV Object
cat("  6. Creating InferCNV object...
")
infercnv_obj <- tryCatch({
     CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                          annotations_file = annotations_file,
                          gene_order_file = gene_order_file_path,
                          ref_group_names = reference_groups)
}, error = function(e) {
    cat(paste("  Error creating InferCNV object:", e$message, "
"))
    return(NULL)
})

if (is.null(infercnv_obj)) {
    cat(paste("  Skipping InferCNV run for stage", stage, "due to object creation error.
"))
    quit(status = 1) # Exit with error status
}

# 7. Run InferCNV
cat("  7. Running InferCNV analysis (this may take a long time)...
")
infercnv_output <- tryCatch({
     infercnv::run(infercnv_obj,
                   cutoff = infercnv_cutoff,
                   out_dir = stage_output_dir,
                   cluster_by_groups = TRUE, # Cluster cells within annotation groups
                   denoise = infercnv_denoise,
                   HMM = infercnv_HMM,
                   num_threads = num_threads,
                   window_length = infercnv_window,
                   analysis_mode = 'samples', # Default
                   output_format = "pdf", # or png
                   write_expr_matrix = TRUE, # Save intermediate matrices
                   plot_steps = FALSE # Don't plot intermediate steps to save time/space
                  )
}, error = function(e) {
    cat(paste("  Error during InferCNV run:", e$message, "
"))
    return(NULL)
})

if (is.null(infercnv_output)) {
    cat(paste("  InferCNV run failed for stage", stage, ". Check logs in:", stage_output_dir, "
"))
    quit(status = 1) # Exit with error status
} else {
    cat(paste("  InferCNV analysis complete for stage", stage, ". Results in:", stage_output_dir, "
"))
}

cat(paste("
======== InferCNV Analysis Script Finished for Stage:", stage, "========
")) 