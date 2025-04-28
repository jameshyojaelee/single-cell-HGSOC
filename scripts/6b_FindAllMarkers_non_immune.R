# ===============================================================================
# FindAllMarkers Analysis and Heatmap for Non-Immune Cells Across Stages
# ===============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(readr)
library(viridis) # For nice color palettes
library(tibble)

# Increase the memory limit for future (important for large datasets)
# Adjust this value based on your system's memory
options(future.globals.maxSize = 50 * 1024^3) # 50GB limit, adjust as needed

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/seurat")
output_dir <- file.path(base_dir, "output/diff_expr/non_immune_all_stages")
plot_dir <- file.path(output_dir, "plots")

# Create output directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Define the path to the *annotated* integrated Seurat object
integrated_object_path <- file.path(input_dir, "integrated_seurat_annotated.rds")

# Check if the annotated integrated object file exists
if (!file.exists(integrated_object_path)) {
  stop(paste("Annotated Seurat object not found at:", integrated_object_path,
             "
Please run script 4_celltype_annotation.R first."))
}

# Load the annotated integrated Seurat object
cat(paste("Loading *annotated* integrated Seurat object from:", integrated_object_path, "
"))
integrated_seurat <- readRDS(integrated_object_path)
cat("Annotated object loaded successfully.
")

# Check for required metadata columns
required_cols <- c("sctype_classification", "tumor_stage")
missing_cols <- required_cols[!required_cols %in% colnames(integrated_seurat@meta.data)]
if (length(missing_cols) > 0) {
  stop(paste("Missing required metadata columns:", paste(missing_cols, collapse=", ")))
}

# Define the target cell type for subsetting
target_cell_type <- "Non-immune cells"

# Check if the target cell type exists
if (!target_cell_type %in% unique(integrated_seurat@meta.data$sctype_classification)) {
    stop(paste("Target cell type '", target_cell_type, "' not found in sctype_classification."))
}

# Subset the Seurat object to keep only the target non-immune cells
cat(paste("Subsetting object to keep only '", target_cell_type, "'...
"))
subset_seurat <- subset(integrated_seurat, subset = sctype_classification == target_cell_type)
cat(paste("Subset created. Number of cells:", ncol(subset_seurat), "
"))

# Set DefaultAssay to RNA before scaling
DefaultAssay(subset_seurat) <- "RNA"

# Check if there are enough cells and stages
if (ncol(subset_seurat) < 10) {
  stop("Too few cells after subsetting for FindAllMarkers analysis.")
}
tumor_stages <- levels(factor(subset_seurat@meta.data$tumor_stage)) # Get all stages present
if (length(tumor_stages) < 2) {
  stop("Insufficient number of stages in the subset for comparison.")
}

# Set cell identities to tumor_stage for FindAllMarkers
cat("Setting cell identities to 'tumor_stage'...
")
Idents(subset_seurat) <- "tumor_stage"

# Ensure the factor levels are set correctly for plotting order if desired
# Example: Order them Normal, IC1, IC2, ...
# stage_order <- c("Normal", "IC1", "IC2", "IC3", "IC4") # Define your desired order
# subset_seurat$tumor_stage <- factor(subset_seurat$tumor_stage, levels = stage_order)
# Idents(subset_seurat) <- "tumor_stage" # Re-set idents after factor ordering

# Set up parallel processing
plan("multicore", workers = 16) 

# --- Run FindAllMarkers ---
cat("Running FindAllMarkers to compare all stages...
")

# Run FindAllMarkers - finds markers for each identity class compared to all remaining cells
# Note: This compares each stage against ALL other stages combined by default
all_stage_markers <- FindAllMarkers(
  subset_seurat,
  assay = "RNA",           # Use RNA assay for biological interpretation
  only.pos = TRUE,         # Keep only positive markers (upregulated in the group)
  test.use = "wilcox",     # Wilcoxon Rank Sum test
  min.pct = 0.1,           # Detect genes expressed in at least 10% of cells in the group
  logfc.threshold = 0.25,  # Minimum log2 fold change
  verbose = TRUE
)

# Check if markers were found
if (nrow(all_stage_markers) == 0) {
  stop("No markers found with the specified thresholds.")
} else {
  cat(paste("Found", nrow(all_stage_markers), "total markers across all stages.
"))
}

# Save the full marker list
output_markers_file <- file.path(output_dir, "FindAllMarkers_non_immune_stages_wilcox.csv")
cat(paste("Saving full marker list to:", output_markers_file, "
"))
write.csv(all_stage_markers, file = output_markers_file, row.names = FALSE)

# --- Select Top Markers for Heatmap ---
cat("Selecting top markers for heatmap...
")
top_n_markers <- 15 # Number of top markers per stage to display
all_stage_markers %>%
  group_by(cluster) %>% # 'cluster' column holds the identity (tumor_stage)
  slice_max(n = top_n_markers, order_by = avg_log2FC) -> top_markers

# Check if top markers were selected
if (nrow(top_markers) == 0) {
  cat("Warning: Could not select top markers. The marker list might be empty or filtering failed.
")
  # Proceed with all markers if top selection fails, or handle appropriately
  top_genes_for_heatmap <- unique(all_stage_markers$gene)
} else {
   cat(paste("Selected top", top_n_markers, "markers per stage based on avg_log2FC.
"))
   top_genes_for_heatmap <- unique(top_markers$gene)
}

# --- Generate Heatmap ---
cat("Generating heatmap using scaled RNA data...
")

# Ensure there are genes to plot
if(length(top_genes_for_heatmap) == 0) {
    stop("No genes selected for heatmap generation.")
}

# Scale data for the selected features using the RNA assay
cat("Scaling data for heatmap features...
")
subset_seurat <- ScaleData(subset_seurat, features = top_genes_for_heatmap, assay = "RNA")

# Create the heatmap using scaled data
heatmap_plot <- DoHeatmap(
  subset_seurat,
  features = top_genes_for_heatmap,
  group.by = "tumor_stage", # Color cells by tumor stage
  assay = "RNA",           # Specify RNA assay
  slot = "scale.data",     # Use scaled data
  label = TRUE,            # Show labels for stages
  size = 5,                # Increase label size
  angle = 0,              # Angle for labels if needed
  hjust = 0.5,               # Center labels horizontally
  # Use a diverging color palette for scaled data
  # group.colors = viridis::viridis(length(tumor_stages)), # Viridis might not be ideal for scaled data
  draw.lines = TRUE        # Draw lines separating groups
) +
# Use a diverging color scale suitable for scaled data (e.g., blue-white-red)
scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
ggtitle("Top Differentially Expressed Genes Across Tumor Stages (Non-Immune Cells)") +
theme(
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
  axis.text.y = element_text(size = 10), # Increase gene label size
  legend.text = element_text(size = 12),      # <-- Increase legend label text size here
  legend.title = element_text(size = 14),     # <-- (Optional) increase legend title size
  legend.position = "right"
)

# Save the heatmap
heatmap_file_path <- file.path(plot_dir, "Heatmap_FindAllMarkers_non_immune_stages_top10.png")
cat(paste("Saving heatmap to:", heatmap_file_path, "
"))
# Increase dimensions and resolution for better quality
ggsave(heatmap_file_path, plot = heatmap_plot, width = 14, height = 10, dpi = 300)

# Reset plan
plan("sequential")

cat("
FindAllMarkers analysis and heatmap generation complete.
")
cat(paste("Full marker list saved in:", output_dir, "
"))
cat(paste("Heatmap saved in:", plot_dir, "
")) 