# ===============================
# Cell Type Annotation with scType
# ===============================

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(HGNChelper)
library(openxlsx)
library(readr)
library(future)
library(ggraph)
library(igraph)
library(tidyverse)
library(data.tree)

# Increase the memory limit for future (important for large datasets)
# Adjust this value based on your system's memory
options(future.globals.maxSize = 200 * 1024^3) # 200GB limit, adjust as needed

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/seurat")
output_dir <- file.path(base_dir, "output/seurat")

# Define the path to the integrated Seurat object
integrated_object_path <- file.path(input_dir, "integrated_seurat.rds")
if (!file.exists(integrated_object_path)) {
  stop(paste("Integrated Seurat object not found at:", integrated_object_path))
}
integrated_seurat <- readRDS(integrated_object_path)

# Load gene set preparation and annotation functions from scType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Set scType DB and tissue type
scType_db <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Immune system" # e.g. Immune system, Pancreas, Liver, etc.

# Prepare gene sets
cat("Preparing gene sets for scType...\n")
gs_list <- gene_sets_prepare(scType_db, tissue)

# Run scType annotation
cat("Running scType annotation...\n")
source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R")
integrated_seurat <- run_sctype(
  integrated_seurat,
  known_tissue_type = tissue,
  custom_marker_file = scType_db,
  name = "sctype_classification",
  plot = FALSE
)

# Check if sctype_classification column exists
if (!"sctype_classification" %in% colnames(integrated_seurat@meta.data)) {
  stop("sctype_classification column not found in Seurat object metadata after run_sctype.")
}

# Check if tumor_stage column exists
if (!"tumor_stage" %in% colnames(integrated_seurat@meta.data)) {
  stop("tumor_stage column not found in Seurat object metadata.")
}

# Reorder tumor_stage levels
integrated_seurat@meta.data$tumor_stage <- factor(
  integrated_seurat@meta.data$tumor_stage,
  levels = c("Normal", "IC2", "IIB", "IIIB", "IIIC")
)

# Fine-tune cell type annotation: merge certain types into 'Non-immune cells'
nonimmune_types <- c(
  "Platelets",
  "Cancer cells",
  "Erythroid-like and erythroid precursor cells",
  "Progenitor cells"
)
integrated_seurat@meta.data$sctype_classification[
  integrated_seurat@meta.data$sctype_classification %in% nonimmune_types
] <- "Non-immune cells"

cat("Updated cell type annotation for fine-tuning. New cell type counts:\n")
print(table(integrated_seurat@meta.data$sctype_classification))

# Export UMAP with cell type annotation (scType)
cat("Exporting UMAP plots...\n")
pdf(file.path(output_dir, "umap_sctype_annotation.pdf"), width = 10, height = 6)
print(
  DimPlot(
    integrated_seurat,
    reduction = "umap",
    group.by = "sctype_classification",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle("UMAP: scType Cell Type Annotation")
)
dev.off()

# Export split UMAP by tumor_stage, annotated by scType
pdf(file.path(output_dir, "umap_sctype_by_tumor_stage.pdf"), width = 28, height = 7)
print(
  DimPlot(
    integrated_seurat,
    reduction = "umap",
    group.by = "sctype_classification",
    split.by = "tumor_stage",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle("UMAP: scType Cell Type Annotation by Tumor Stage")
)
dev.off()

# Calculate cell type counts and percentages by tumor stage
cat("\nCalculating cell type distributions...\n")
cell_counts <- integrated_seurat@meta.data %>%
  group_by(tumor_stage, sctype_classification) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  ungroup()

# Calculate total cells per tumor stage for percentage calculation
stage_totals <- integrated_seurat@meta.data %>%
  group_by(tumor_stage) %>%
  summarise(stage_total = n(), .groups = 'drop')

# Join counts with totals and calculate percentages
cell_distribution <- cell_counts %>%
  left_join(stage_totals, by = "tumor_stage") %>%
  mutate(
    percentage = round(cell_count / stage_total * 100, 2),
    percentage_label = sprintf("%.2f%%", percentage)
  ) %>%
  arrange(
    factor(tumor_stage, levels = c("Normal", "IC2", "IIB", "IIIB", "IIIC")),
    desc(cell_count)
  )

# Create a formatted table with counts and percentages
distribution_table <- cell_distribution %>%
  select(
    `Tumor Stage` = tumor_stage,
    `Cell Type` = sctype_classification,
    `Cell Count` = cell_count,
    `Percentage` = percentage_label,
    `Stage Total` = stage_total
  )

# Export results
cat("Exporting cell type distribution results...\n")
write.csv(
  distribution_table,
  file = file.path(output_dir, "celltype_distribution_by_stage.csv"),
  row.names = FALSE
)

# Print summary to console
cat("\nCell type distribution summary:\n")
print(distribution_table, n = Inf)

# Create a summary of total cells per stage
stage_summary <- stage_totals %>%
  arrange(factor(tumor_stage, levels = c("Normal", "IC2", "IIB", "IIIB", "IIIC"))) %>%
  rename(
    `Tumor Stage` = tumor_stage,
    `Total Cells` = stage_total
  )

cat("\nTotal cells per tumor stage:\n")
print(stage_summary, n = Inf)

cat("\nResults exported to:", file.path(output_dir, "celltype_distribution_by_stage.csv"), "\n")

# Save the updated Seurat object with annotations
annotated_object_path <- file.path(output_dir, "integrated_seurat_annotated.rds")
cat(paste("\nSaving annotated Seurat object to:", annotated_object_path, "\n"))
saveRDS(integrated_seurat, file = annotated_object_path)

cat("Script complete!\n")
