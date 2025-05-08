# ======================================================================
# Differential Expression Analysis for T Cells Across Stages
# ======================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(future)
library(readr)

# Increase the memory limit for future (important for large datasets)
# Adjust this value based on your system's memory
options(future.globals.maxSize = 50 * 1024^3) # 50GB limit, adjust as needed

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/seurat")
output_dir <- file.path(base_dir, "output/diff_expr/t_cells")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the path to the *annotated* integrated Seurat object
integrated_object_path <- file.path(input_dir, "integrated_seurat_annotated.rds")

# Check if the annotated integrated object file exists
if (!file.exists(integrated_object_path)) {
  stop(paste("Annotated Seurat object not found at:", integrated_object_path, 
             "\nPlease run script 4_celltype_annotation.R first."))
}

# Load the annotated integrated Seurat object
cat(paste("Loading *annotated* integrated Seurat object from:", integrated_object_path, "\n"))
integrated_seurat <- readRDS(integrated_object_path)
cat("Annotated object loaded successfully.\n")

# Check for required metadata columns
required_cols <- c("sctype_classification", "tumor_stage")
missing_cols <- required_cols[!required_cols %in% colnames(integrated_seurat@meta.data)]
if (length(missing_cols) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing_cols, collapse=", ")))
}

# Define the target cell types for subsetting (all T cell subtypes)
t_cell_types <- c("Naive CD4+ T cells", "Memory CD4+ T cells", "CD4+ NKT-like cells", "CD8+ NKT-like cells")

# Check if any of the target cell types exist
present_t_cell_types <- intersect(t_cell_types, unique(integrated_seurat@meta.data$sctype_classification))
if (length(present_t_cell_types) == 0) {
    stop(paste("None of the T cell subtypes found in sctype_classification. Checked:", paste(t_cell_types, collapse=", ")))
}

# Subset the Seurat object to keep only the target T cell subtypes
cat(paste("Subsetting object to keep only T cell subtypes:", paste(present_t_cell_types, collapse=", "), "...\n"))
subset_seurat <- subset(integrated_seurat, subset = sctype_classification %in% present_t_cell_types)
cat(paste("Subset created. Number of cells:", ncol(subset_seurat), "\n"))

# Check if there are enough cells and stages for DEG
if (ncol(subset_seurat) < 10) {
  stop("Too few cells after subsetting for DEG analysis.")
}
tumor_stages <- unique(subset_seurat@meta.data$tumor_stage)
if (!"Normal" %in% tumor_stages || length(tumor_stages) < 2) {
  stop("'Normal' stage not found or insufficient number of stages in the subset for comparison.")
}

# Set cell identities to tumor_stage for DEG analysis
cat("Setting cell identities to 'tumor_stage' for DEG...\n")
Idents(subset_seurat) <- "tumor_stage"

# Get stages to compare against Normal (excluding Normal itself)
comparison_stages <- setdiff(tumor_stages, "Normal")

# Set up parallel processing
plan("multicore", workers = 8) # Adjust workers based on your system

# Perform DEG analysis for each stage vs. Normal
cat("Starting DEG analysis (comparing each stage vs. Normal)...")

all_markers <- list()

for (stage in comparison_stages) {
  cat(paste("\nComparing", stage, "vs. Normal...\n"))
  
  # Check if both groups have enough cells
  cells_stage <- WhichCells(subset_seurat, idents = stage)
  cells_normal <- WhichCells(subset_seurat, idents = "Normal")
  
  if (length(cells_stage) < 3 || length(cells_normal) < 3) {
      cat(paste("Skipping comparison for", stage, "vs Normal due to insufficient cells in one or both groups.\n"))
      next
  }
  
  tryCatch({
    markers <- FindMarkers(
      subset_seurat,
      ident.1 = stage,
      ident.2 = "Normal",
      assay = "RNA", # Explicitly use the RNA assay for DEG
      test.use = "wilcox", # Wilcoxon Rank Sum test
      logfc.threshold = 0.25, # Consider adjusting if needed
      min.pct = 0.1,        # Consider adjusting if needed
      verbose = TRUE
    )
    
    if (nrow(markers) > 0) {
      cat(paste("Found", nrow(markers), "markers for", stage, "vs. Normal.\n"))
      markers$gene <- rownames(markers)
      markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      all_markers[[paste0(stage, "_vs_Normal")]] <- markers
      
      # Save results for this comparison
      output_file_path <- file.path(output_dir, paste0("DEG_", stage, "_vs_Normal_wilcox.csv"))
      cat(paste("Saving results to:", output_file_path, "\n"))
      write.csv(markers, file = output_file_path, row.names = FALSE)
      
    } else {
      cat(paste("No significant markers found for", stage, "vs. Normal with current thresholds.\n"))
    }
    
  }, error = function(e) {
    cat(paste("Error during FindMarkers for", stage, "vs. Normal:", e$message, "\n"))
  })
}

# Reset plan
plan("sequential")

cat("\nDEG analysis complete for T cells.\n")
cat(paste("Results saved in:", output_dir, "\n")) 