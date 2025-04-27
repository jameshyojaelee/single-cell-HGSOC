# ======================================================
# Cell-Cell Interaction Analysis using CellChat
# ======================================================

# --- Prerequisites ---
# 1. Install CellChat (see vignette or GitHub)
#    devtools::install_github("sqjin/CellChat")

# Load required libraries
library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(future)

# Increase memory limit if needed (CellChat can be memory intensive)
options(future.globals.maxSize = 200 * 1024^3) # 80GB limit, adjust as needed

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/seurat") # Directory where annotated object is saved
output_dir <- file.path(base_dir, "output/cellchat")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Input file
annotated_object_path <- file.path(input_dir, "integrated_seurat_annotated.rds")

# Check if the annotated object file exists
if (!file.exists(annotated_object_path)) {
  stop(paste("Annotated Seurat object not found at:", annotated_object_path,
             "\nPlease run script 4_celltype_annotation.R first."))
}

# Parameters
cell_identity_col <- "sctype_classification" # Metadata column for cell type labels
condition_col <- "tumor_stage" # Metadata column for conditions/stages
min_cells_per_type <- 10 # Minimum cells required for a cell type to be included
min_cells_per_dataset <- 50 # Minimum cells required for a stage/dataset

# --- Load Data and Database ---

cat(paste("Loading annotated Seurat object from:", annotated_object_path, "\n"))
seurat_obj <- readRDS(annotated_object_path)
cat("Annotated object loaded successfully.\n")

# Check required metadata
if (!cell_identity_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required cell identity column not found:", cell_identity_col))
}
if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required condition column not found:", condition_col))
}

# Load CellChat database for Human
cat("Loading CellChatDB...")
cellchat_db <- CellChatDB.human
cat(" Done.\n")

# --- Analysis Loop for Each Stage ---

stages <- levels(seurat_obj@meta.data[[condition_col]])
cellchat_list <- list()
processed_stages <- c()

cat(paste("Starting CellChat analysis for each stage:", paste(stages, collapse=", "), "\n"))

for (stage in stages) {
  cat(paste("\nProcessing stage:", stage, "\n"))
  
  # Subset Seurat object for the current stage
  seurat_subset <- subset(seurat_obj, subset = !!sym(condition_col) == stage)
  
  if (ncol(seurat_subset) < min_cells_per_dataset) {
    cat(paste("  Skipping stage", stage, ": Less than", min_cells_per_dataset, "cells total.\n"))
    next
  }
  
  # Prepare input data for CellChat
  data_input <- GetAssayData(seurat_subset, assay = "RNA", slot = "data") # Normalized data
  meta <- seurat_subset@meta.data
  cell_labels <- meta[[cell_identity_col]]
  
  # Filter cell types with insufficient cells
  cell_counts <- table(cell_labels)
  valid_cell_types <- names(cell_counts[cell_counts >= min_cells_per_type])
  if (length(valid_cell_types) < 2) {
     cat(paste("  Skipping stage", stage, ": Less than 2 cell types with >=", min_cells_per_type, "cells.\n"))
     next
  }
  keep_cells <- rownames(meta[meta[[cell_identity_col]] %in% valid_cell_types, ])
  data_input <- data_input[, keep_cells]
  meta <- meta[keep_cells, ]
  cat(paste("  Kept", length(valid_cell_types), "cell types with >=", min_cells_per_type, "cells.\n"))

  # Create CellChat object
  tryCatch({
    cellchat <- createCellChat(object = data_input, meta = meta, group.by = cell_identity_col)
    cellchat <- addMeta(cellchat, meta = meta, meta.name = "metadata")
    cellchat <- setIdent(cellchat, ident.use = cell_identity_col)
    groupSize <- as.numeric(table(cellchat@idents)) # Number of cells in each type
    
    # Set ligand-receptor database
    cellchat@DB <- cellchat_db
    
    # Preprocessing
    cat("  Preprocessing data (identifyOverExpressedGenes, projectData)...")
    cellchat <- subsetData(cellchat) # Subset the expression data for speed
    plan("multicore", workers = 4) # Use future for parallelization
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    plan("sequential")
    cat(" Done.\n")
    
    # Infer communication network
    cat("  Inferring communication network (computeCommunProb, aggregateNet)...")
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE) # Use raw data
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cat(" Done.\n")
    
    # --- Save Object and Generate Plots for this stage ---
    stage_output_dir <- file.path(output_dir, stage)
    if (!dir.exists(stage_output_dir)) dir.create(stage_output_dir)
    
    # Save CellChat object
    saveRDS(cellchat, file = file.path(stage_output_dir, paste0("cellchat_", stage, ".rds")))
    cat(paste("  Saved CellChat object for stage", stage, "\n"))
    
    # Create plots
    pdf(file.path(stage_output_dir, paste0("plots_", stage, ".pdf")), width = 10, height = 8)
    
    # Interaction weights/strength heatmap
    print(netAnalysis_heatmap(cellchat, measure = "weight", 
                              title = paste("Interaction Weights -", stage)))
                              
    # Circle plot
    print(netVisual_circle(cellchat@net$weight, 
                           vertex.weight = groupSize, weight.scale = T, 
                           label.edge= F, title.name = paste("Interaction Weights/Number -", stage)))
                           
    # Add more plots as needed (e.g., pathway level, specific L-R pairs)
    # netVisual_bubble(...), plotGeneExpression(...), etc.
    
    dev.off()
    cat(paste("  Saved plots for stage", stage, "\n"))
    
    # Store the object for potential comparison later
    cellchat_list[[stage]] <- cellchat
    processed_stages <- c(processed_stages, stage)
    
  }, error = function(e) {
    cat(paste("  Error processing stage", stage, ":", e$message, "\n"))
  })
}

# --- Comparative Analysis (Optional) ---

if (length(cellchat_list) >= 2) {
  cat("\nPerforming comparative analysis across stages...\n")
  comparison_output_dir <- file.path(output_dir, "comparison")
  if (!dir.exists(comparison_output_dir)) dir.create(comparison_output_dir)

  tryCatch({
    # Merge objects
    cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    saveRDS(cellchat_merged, file = file.path(comparison_output_dir, "cellchat_merged.rds"))
    
    # --- Generate Comparison Plots ---
    pdf(file.path(comparison_output_dir, "plots_comparison.pdf"), width = 12, height = 10)
    
    # Compare total interactions
    print(compareInteractions(cellchat_merged, show.legend = F, group = names(cellchat_list)))
    
    # Compare interactions heatmap
    print(netAnalysis_heatmap(cellchat_merged, measure = "weight", 
                              title = "Comparison: Interaction Weights"))
                              
    # Differential interactions / Network centrality
    # This requires more specific choices, refer to CellChat comparison vignette
    # Example: Rank networks based on centrality
    # num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)}) ...
    # netAnalysis_rankCentralityComparison(...) 
    
    dev.off()
    cat("  Saved comparison plots.\n")

  }, error = function(e) {
      cat(paste("  Error during comparative analysis:", e$message, "\n"))
  })
} else {
  cat("\nSkipping comparative analysis: Need at least two successfully processed stages.\n")
}

cat("\nCellChat analysis complete.\n")
cat(paste("Results saved in:", output_dir, "\n")) 