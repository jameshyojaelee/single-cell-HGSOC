# ======================================================
# Cell-Cell Interaction Analysis using CellChat
# ======================================================

# Load required libraries
library(CellChat, lib="/gpfs/commons/home/jameslee/R/x86_64-pc-linux-gnu-library/4.4/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(future)
library(stringr) # Added for str_remove

# Increase memory limit if needed (CellChat can be memory intensive)
options(future.globals.maxSize = 200 * 1024^3)

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
original_cell_identity_col <- "sctype_classification" # Original fine-grained cell type labels
grouped_cell_identity_col <- "cell_group"         # New column for broader cell groups
cell_identity_col <- grouped_cell_identity_col      # *** Use grouped identities for CellChat ***
condition_col <- "tumor_stage" # Metadata column for conditions/stages
min_cells_per_type <- 10 # Minimum cells required for a *grouped* cell type to be included
min_cells_per_dataset <- 50 # Minimum cells required for a stage/dataset

# --- Load Data and Database ---

cat(paste("Loading annotated Seurat object from:", annotated_object_path, "\n"))
seurat_obj <- readRDS(annotated_object_path)
cat("Annotated object loaded successfully.\n")

# Check required metadata
if (!original_cell_identity_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required original cell identity column not found:", original_cell_identity_col))
}
if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Required condition column not found:", condition_col))
}

# --- Create Broader Cell Groups ---
cat("Creating broader cell groups...\n")
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(
    !!sym(grouped_cell_identity_col) := case_when(
      !!sym(original_cell_identity_col) %in% c("Naive CD4+ T cells", "Memory CD4+ T cells",
                                               "CD8+ NKT-like cells", "CD4+ NKT-like cells") ~ "T Cells",
      !!sym(original_cell_identity_col) %in% c("Naive B cells", "Pre-B cells") ~ "B Cells",
      !!sym(original_cell_identity_col) %in% c("Macrophages", "Non-classical monocytes",
                                               "Myeloid Dendritic cells") ~ "Myeloid",
      !!sym(original_cell_identity_col) == "Neutrophils" ~ "Neutrophils",
      !!sym(original_cell_identity_col) == "ISG expressing immune cells" ~ "Other Immune",
      !!sym(original_cell_identity_col) == "Endothelial" ~ "Endothelial",
      !!sym(original_cell_identity_col) == "Non-immune cells" ~ "Other Non-immune",
      !!sym(original_cell_identity_col) == "Unknown" ~ "Unknown",
      TRUE ~ "Other" # Catch any unexpected types
    )
  )
# Convert the new group column to a factor for consistency
seurat_obj@meta.data[[grouped_cell_identity_col]] <- factor(seurat_obj@meta.data[[grouped_cell_identity_col]])
cat(paste("  Added metadata column:", grouped_cell_identity_col, "\n"))
print(paste("  Grouped cell types created:"))
print(table(seurat_obj@meta.data[[grouped_cell_identity_col]]))

# Load CellChat database for Human
cat("Loading CellChatDB...")
cellchat_db <- CellChatDB.human
ppi.human <- PPI.human
cat(" Done.\n")

# --- Analysis Loop for Each Stage ---

# Ensure stage levels are correctly ordered if necessary (e.g., based on script 5)
stage_levels_ordered <- c("Normal", "IC2", "IIB", "IIIB", "IIIC")
seurat_obj@meta.data[[condition_col]] <- factor(seurat_obj@meta.data[[condition_col]], levels = stage_levels_ordered)
seurat_obj <- subset(seurat_obj, subset = !is.na(!!sym(condition_col))) # Remove cells with NA stage

stages <- levels(seurat_obj@meta.data[[condition_col]])
cellchat_list <- list()
processed_stages <- c()

cat(paste("Starting CellChat analysis using grouped cell types for each stage:", paste(stages, collapse=", "), "\n"))

for (stage in stages) {
  cat(paste("\nProcessing stage:", stage, "\n"))
  
  # Subset Seurat object for the current stage
  seurat_subset <- subset(seurat_obj, subset = !!sym(condition_col) == stage)
  
  if (ncol(seurat_subset) < min_cells_per_dataset) {
    cat(paste("  Skipping stage", stage, ": Less than", min_cells_per_dataset, "cells total.\n"))
    next
  }
  
  # Prepare input data for CellChat
  data_input <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
  meta <- seurat_subset@meta.data
  cell_labels <- meta[[cell_identity_col]]
  
  # Filter cell types with insufficient cells
  cell_counts <- table(cell_labels)
  valid_cell_types <- names(cell_counts[cell_counts >= min_cells_per_type])
  if (length(valid_cell_types) < 2) {
     cat(paste("  Skipping stage", stage, ": Less than 2 *grouped* cell types with >=", min_cells_per_type, "cells.\n"))
     next
  }
  keep_cells <- rownames(meta[meta[[cell_identity_col]] %in% valid_cell_types, ])
  data_input_filtered <- data_input[, keep_cells]
  meta_filtered <- meta[keep_cells, ]
  cat(paste("  Kept", length(valid_cell_types), "*grouped* cell types with >=", min_cells_per_type, "cells:", paste(valid_cell_types, collapse=", "), "\n"))

  # Ensure the identity column in the filtered metadata is a factor with correct levels
  meta_filtered[[cell_identity_col]] <- factor(meta_filtered[[cell_identity_col]], levels = valid_cell_types)

  # Create CellChat object
  tryCatch({
    cellchat <- createCellChat(object = data_input_filtered, meta = meta_filtered, group.by = cell_identity_col)
    # cellchat <- addMeta(cellchat, meta = meta_filtered, meta.name = "metadata") # No longer needed as group.by handles it? Check CellChat docs. Let's keep addMeta for safety.
    # cellchat <- addMeta(cellchat, meta = meta_filtered, meta.name = "metadata") # Removed: Likely redundant with group.by
    # cellchat <- setIdent(cellchat, ident.use = cell_identity_col) # Removed: createCellChat with group.by should handle this
    
    # Check idents *immediately* after creation with group.by
    cat("    DEBUG: Checking idents immediately after createCellChat(group.by = ...)...\n")
    if (is.null(cellchat@idents) || length(levels(cellchat@idents)) == 0) {
        stop("Error: cellchat@idents are NULL or have zero levels after createCellChat with group.by.")
    } else {
        cat(paste("      Levels set correctly in cellchat@idents by createCellChat:", paste(levels(cellchat@idents), collapse=", "), "\n"))
    }
    
    # Recalculate groupSize based on the idents set by createCellChat
    groupSize <- as.numeric(table(cellchat@idents)) # Number of cells in each grouped type

    # <<< --- ADD IDENTS CHECK HERE --- >>>
    # This check is now redundant with the one above, but keep for sanity
    cat("    DEBUG: Re-checking idents after groupSize calculation (using grouped types)...\n")
    if (is.null(cellchat@idents) || length(levels(cellchat@idents)) == 0) {
        # This error should theoretically not be reached if the first check passed
        stop("Error: cellchat@idents became NULL or lost levels after groupSize calculation.") 
    } else if (length(levels(cellchat@idents)) != length(valid_cell_types)) {
        cat(paste("      Warning: Number of levels in idents (", length(levels(cellchat@idents)), ") does not match number of valid cell types (", length(valid_cell_types), ").\n"))
        cat(paste("      Levels in idents:", paste(levels(cellchat@idents), collapse=", "), "\n"))
        cat(paste("      Valid cell types defined earlier:", paste(valid_cell_types, collapse=", "), "\n"))
    } else {
        cat(paste("      Levels still set correctly in cellchat@idents:", paste(levels(cellchat@idents), collapse=", "), "\n"))
    }
    cat("    DEBUG: Idents check complete.\n")
    # <<< --- END IDENTS CHECK --- >>>
    
    # Set ligand-receptor database
    cellchat@DB <- cellchat_db
    
    # Preprocessing
    cat("  Preprocessing data (subsetData)...")
    cellchat <- subsetData(cellchat) # Subset the expression data for speed
    cat(" Done subsetData.\n")

    # <<< --- ADD DEBUG CHECKS HERE --- >>>
    cat("    DEBUG: Checking genes before identifyOverExpressedGenes...\n")
    # Check original input genes against DB
    input_genes <- rownames(data_input) 
    db_genes_ligand <- unique(cellchat@DB$interaction$ligand)
    db_genes_receptor <- unique(unlist(strsplit(cellchat@DB$interaction$receptor, "&")))
    db_genes_interaction <- unique(c(db_genes_ligand, db_genes_receptor))
    common_genes_input_db <- intersect(input_genes, db_genes_interaction)
    cat(paste("      Genes in input data:", length(input_genes), "\n"))
    cat(paste("      Genes in DB interactions:", length(db_genes_interaction), "\n"))
    cat(paste("      Common genes (Input vs DB interactions):", length(common_genes_input_db), "\n"))
    if (length(common_genes_input_db) < 50) { # Arbitrary low threshold
        cat("      WARNING: Very few common genes found between input and DB!\n")
    }
    
    # Check genes remaining after subsetData
    genes_after_subset <- rownames(cellchat@data.signaling)
    cat(paste("      Genes remaining in cellchat@data.signaling (after subsetData):", length(genes_after_subset), "\n"))
    if (length(genes_after_subset) == 0) {
        stop("Error: No genes remaining after subsetData. Cannot proceed.")
    }
    cat(paste("      Dimensions of cellchat@data.signaling:", paste(dim(cellchat@data.signaling), collapse=" x "), "\n"))
    cat("    DEBUG: Checks complete.\n")
    # <<< --- END DEBUG CHECKS --- >>>
    
    cat("  Preprocessing data (identifyOverExpressedGenes, projectData)...")
    plan("multicore", workers = 4) # Use future for parallelization
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    plan("sequential")
    cat(" Done identify/project.\n")
    
    # Infer communication network
    cat("  Inferring communication network (computeCommunProb, aggregateNet)...")
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE) # Use raw data
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cat(" Done.\n")
    
    # --- Save Object for this stage ---
    stage_output_dir <- file.path(output_dir, stage)
    if (!dir.exists(stage_output_dir)) dir.create(stage_output_dir)
    
    # Save CellChat object
    saveRDS(cellchat, file = file.path(stage_output_dir, paste0("cellchat_grouped_", stage, ".rds")))
    cat(paste("  Saved CellChat object (grouped) for stage", stage, "\n"))
    
    # Store the object for potential comparison later
    cellchat_list[[stage]] <- cellchat
    processed_stages <- c(processed_stages, stage)
    
  }, error = function(e) {
    cat(paste("  Error processing stage", stage, ":", e$message, "\n"))
  })
}
########################################################


# --- Comparative Analysis (Optional) ---

if (length(cellchat_list) >= 2) {
  cat("\nPerforming comparative analysis across stages (using grouped types)...\n")
  comparison_output_dir <- file.path(output_dir, "comparison_grouped")
  if (!dir.exists(comparison_output_dir)) dir.create(comparison_output_dir)

  tryCatch({
    # Merge objects
    cat("  Merging CellChat objects...")
    cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    saveRDS(cellchat_merged, file = file.path(comparison_output_dir, "cellchat_merged_grouped.rds"))
    cat(" Merged and saved.\n")

  }, error = function(e) {
      cat(paste("  Error during comparative analysis (grouped):", e$message, "\n"))
  })
} else {
  cat("\nSkipping comparative analysis: Need at least two successfully processed stages.\n")
}

cat("\nCellChat analysis complete.\n")
cat(paste("Results (using grouped cell types) saved in:", output_dir, "\nRun script 11_visualize_cellchat.R to generate plots.\n"))