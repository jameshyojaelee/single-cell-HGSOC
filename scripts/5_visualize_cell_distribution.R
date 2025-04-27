# =======================================================
# Visualize Cell Type Distribution Changes Across Stages
# =======================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(tidyr)
library(RColorBrewer)
library(stringr)

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_file <- file.path(base_dir, "output/seurat/celltype_distribution_by_stage.csv")
output_dir <- file.path(base_dir, "output/plots")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check if input file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Load the data
cat(paste("Loading data from:", input_file, "\n"))
cell_dist <- read_csv(input_file, show_col_types = FALSE)

# --- Data Preprocessing ---
cat("Preprocessing data...\n")

# Define the desired order for tumor stages
stage_levels <- c("Normal", "IC2", "IIB", "IIIB", "IIIC")

# Convert Percentage string to numeric and order Tumor Stage
cell_dist <- cell_dist %>%
  mutate(
    Percentage_Num = as.numeric(str_remove(`Percentage`, "%$")),
    `Tumor Stage` = factor(`Tumor Stage`, levels = stage_levels)
  ) %>%
  filter(!is.na(`Tumor Stage`)) # Remove rows with potential NA stages if any

# --- Create Broader Cell Groups ---
cell_dist <- cell_dist %>%
  mutate(
    `Cell Group` = case_when(
      `Cell Type` %in% c("Naive CD4+ T cells", "Memory CD4+ T cells", 
                         "CD8+ NKT-like cells", "CD4+ NKT-like cells") ~ "T Cells",
      `Cell Type` %in% c("Naive B cells", "Pre-B cells") ~ "B Cells",
      `Cell Type` %in% c("Macrophages", "Non-classical monocytes", 
                         "Myeloid Dendritic cells") ~ "Myeloid",
      `Cell Type` == "Neutrophils" ~ "Neutrophils",
      `Cell Type` == "ISG expressing immune cells" ~ "Other Immune",
      `Cell Type` == "Endothelial" ~ "Endothelial",
      `Cell Type` == "Non-immune cells" ~ "Other Non-immune",
      `Cell Type` == "Unknown" ~ "Unknown",
      TRUE ~ "Other" # Catch any unexpected types
    )
  )

# Aggregate percentages by the new Cell Group
cell_dist_grouped <- cell_dist %>%
  group_by(`Tumor Stage`, `Cell Group`) %>%
  summarise(Percentage = sum(Percentage_Num), .groups = 'drop')

# Define color palettes
# Palette for broader groups (adjust if more than 9 groups)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
n_groups <- length(unique(cell_dist_grouped$`Cell Group`))
group_palette <- colorRampPalette(brewer.pal(min(n_groups, 9), "Set1"))(n_groups)

# Palette for fine types (can become large)
n_fine_types <- length(unique(cell_dist$`Cell Type`))
fine_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_fine_types)

# --- Generate Plots ---
cat("Generating plots...\n")

# 1. Stacked Bar Plot (Grouped Cell Types)
plot1 <- ggplot(cell_dist_grouped, aes(x = `Tumor Stage`, y = Percentage, fill = `Cell Group`)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Cell Type Group Proportions Across Tumor Stages",
    x = "Tumor Stage",
    y = "Percentage of Total Cells per Stage",
    fill = "Cell Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "1_stacked_bar_grouped.pdf"), plot1, width = 8, height = 6)
cat("Saved: 1_stacked_bar_grouped.pdf\n")

# 2. Faceted Bar Plot (Grouped Cell Types)
plot2 <- ggplot(cell_dist_grouped, aes(x = `Tumor Stage`, y = Percentage / 100, fill = `Cell Group`)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ `Cell Group`, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Percentage of Each Cell Group Across Tumor Stages",
    x = "Tumor Stage",
    y = "Percentage of Total Cells per Stage"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.background = element_rect(fill="grey90"),
        panel.spacing = unit(1, "lines"))

ggsave(file.path(output_dir, "2_faceted_bar_grouped.pdf"), plot2, width = 10, height = 8)
cat("Saved: 2_faceted_bar_grouped.pdf\n")

# 3. Line Plot (Grouped Cell Types)
plot3 <- ggplot(cell_dist_grouped, aes(x = `Tumor Stage`, y = Percentage / 100, 
                                       group = `Cell Group`, color = `Cell Group`)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = group_palette) +
  labs(
    title = "Trend of Cell Group Percentages Across Tumor Stages",
    x = "Tumor Stage",
    y = "Percentage of Total Cells per Stage",
    color = "Cell Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "3_line_plot_grouped.pdf"), plot3, width = 9, height = 6)
cat("Saved: 3_line_plot_grouped.pdf\n")

# 4. Heatmap (Grouped Cell Types)
plot4 <- ggplot(cell_dist_grouped, aes(x = `Tumor Stage`, y = `Cell Group`, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Percentage") +
  labs(
    title = "Heatmap of Cell Group Percentages Across Tumor Stages",
    x = "Tumor Stage",
    y = "Cell Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "4_heatmap_grouped.pdf"), plot4, width = 8, height = 6)
cat("Saved: 4_heatmap_grouped.pdf\n")

# 5. Stacked Bar Plot (Fine Cell Types)
# Ensure factor levels are set for fine types if specific order is needed
cell_dist$`Cell Type` <- factor(cell_dist$`Cell Type`)

plot5 <- ggplot(cell_dist, aes(x = `Tumor Stage`, y = Percentage_Num, fill = `Cell Type`)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = fine_palette) +
  labs(
    title = "Detailed Cell Type Proportions Across Tumor Stages",
    x = "Tumor Stage",
    y = "Percentage of Total Cells per Stage",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=8))

ggsave(file.path(output_dir, "5_stacked_bar_fine.pdf"), plot5, width = 10, height = 7)
cat("Saved: 5_stacked_bar_fine.pdf\n")

cat("\nAll plots saved to:", output_dir, "\n") 