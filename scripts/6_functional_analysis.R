# =======================================================
# Functional Analysis of Immune Cell States in HGSOC
# =======================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(RColorBrewer)
library(stringr)
library(patchwork) # For combining plots
library(ggpubr)    # For statistical comparisons

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
seurat_object_path <- file.path(base_dir, "output/seurat/integrated_seurat_annotated.rds")
output_dir <- file.path(base_dir, "output/plots/functional_analysis")
metadata_column <- "sctype_classification" # Column containing cell type annotations
stage_column <- "tumor_stage" # Assuming 'Tumor_Stage' based on previous script, confirm if different

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check if Seurat object file exists
if (!file.exists(seurat_object_path)) {
  stop(paste("Seurat object file not found:", seurat_object_path))
}

# Load the Seurat object
cat(paste("Loading Seurat object from:", seurat_object_path, "\n"))
hgsoc_seurat <- readRDS(seurat_object_path)
cat("Seurat object loaded.\n")

# --- Define Cell Types of Interest ---
cat("Defining cell types of interest...\n")
cd4_types <- c("Naive CD4+ T cells", "Memory CD4+ T cells", "CD4+ NKT-like cells")
cd8_types <- c("CD8+ NKT-like cells")
# Note: Using "CD8+ NKT-like cells" for CD8 analysis.
macrophage_types <- c("Macrophages", "Non-classical monocytes")
# Note: No specific 'NK cells' found. Analyzing 'CD8+ NKT-like cells' and 'ISG expressing immune cells' for NK markers.
nk_types <- c("CD8+ NKT-like cells", "ISG expressing immune cells")

# Combine T cell types for easier subsetting
t_cell_types <- c(cd4_types, cd8_types)

# --- Verify Metadata and Genes ---
# Check if metadata columns exist
if (!metadata_column %in% colnames(hgsoc_seurat@meta.data)) {
  stop(paste("Metadata column '", metadata_column, "' not found in Seurat object."))
}
# Check if stage column exists (optional, for grouping plots)
if (!stage_column %in% colnames(hgsoc_seurat@meta.data)) {
    warning(paste("Stage column '", stage_column, "' not found. Plots will not be grouped by stage."))
    stage_column <- NULL # Set to NULL if not found
} else {
    # Ensure stage column is a factor with desired levels if needed (similar to script 5)
    # We define levels and palette later, only if stage_column is confirmed valid for plotting
    # stage_levels <- c("Normal", "IC2", "IIB", "IIIB", "IIIC") # Moved down
    # hgsoc_seurat[[stage_column]] <- factor(hgsoc_seurat[[stage_column, drop=TRUE]], levels = stage_levels)
    pass # Placeholder, just confirm the column exists here
}


# Define genes for analysis (adjust symbols if necessary)
exhaustion_markers <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
nk_markers <- c("NCAM1", "NKG7", "GNLY", "KLRF1", "KLRD1", "FCGR3A") # CD56, NKG7, Granulysin, NKp80, CD94, CD16
m1_genes <- list(M1 = c("TNF", "IL6", "IL1B", "CXCL9", "CXCL10"))
m2_genes <- list(M2 = c("CD163", "MRC1", "IL10", "CCL18", "CCL22"))

# Check if genes exist in the object (using RNA assay)
all_genes_list <- c(exhaustion_markers, nk_markers, unlist(m1_genes), unlist(m2_genes))
present_genes <- all_genes_list[all_genes_list %in% rownames(hgsoc_seurat[["RNA"]])]
missing_genes <- all_genes_list[!all_genes_list %in% rownames(hgsoc_seurat[["RNA"]])]

if (length(missing_genes) > 0) {
    warning("The following genes were not found in the Seurat object (RNA assay) and will be excluded: ", paste(missing_genes, collapse=", "))
}

# Update gene lists to only include present genes
exhaustion_markers <- exhaustion_markers[exhaustion_markers %in% present_genes]
nk_markers <- nk_markers[nk_markers %in% present_genes]
m1_genes[[1]] <- m1_genes[[1]][m1_genes[[1]] %in% present_genes]
m2_genes[[1]] <- m2_genes[[1]][m2_genes[[1]] %in% present_genes]


# --- 1. T Cell Exhaustion Analysis ---
cat("\n--- Starting T cell exhaustion analysis ---\n")

if (length(exhaustion_markers) > 0) {
    # Subset Seurat object for T cells using original types
    t_cell_subset <- subset(hgsoc_seurat, subset = !!sym(metadata_column) %in% t_cell_types)
    
    # Define simplified names
    cd4_target_name <- "CD4+ T cells"
    cd8_target_name <- "CD8+ T cells"
    
    # Add simplified names to metadata - apply directly to the metadata dataframe
    t_cell_subset@meta.data <- t_cell_subset@meta.data %>% 
        mutate(simplified_cell_type = case_when(
            !!sym(metadata_column) %in% cd4_types ~ cd4_target_name,
            !!sym(metadata_column) %in% cd8_types ~ cd8_target_name,
            TRUE ~ as.character(!!sym(metadata_column)) # Should not happen based on subsetting
        ))
    # Ensure it's a factor for plotting
    t_cell_subset$simplified_cell_type <- factor(t_cell_subset$simplified_cell_type)

    # Only generate plots if stage column is available
    if (!is.null(stage_column)) {
        cat("Generating T cell exhaustion plots grouped by stage (simplified names)...\n")

        # Define stage levels and palette here, ensuring they exist in this scope
        stage_levels <- c("Normal", "IC2", "IIB", "IIIB", "IIIC")
        # Ensure the column is factored correctly within the subset *before* plotting
        t_cell_subset[[stage_column]] <- factor(t_cell_subset[[stage_column, drop=TRUE]], levels = stage_levels)
        
        n_stages <- length(stage_levels)
        # Change palette to Yellow-Orange-Red
        stage_palette <- RColorBrewer::brewer.pal(max(3, n_stages), "YlOrRd")[1:n_stages]
        names(stage_palette) <- stage_levels

        # Plot 1: Grouped Violin Plot by Stage
        plot_vln_exhaustion_stage_list <- VlnPlot(t_cell_subset, 
                                                  features = exhaustion_markers, 
                                                  group.by = "simplified_cell_type", # Use simplified names
                                                  split.by = stage_column,
                                                  pt.size = 0, 
                                                  combine = FALSE, 
                                                  log = TRUE,
                                                  raster = FALSE) # Disable rasterization
        # Apply theme modifications and color palette to each plot in the list
        plot_vln_exhaustion_stage <- lapply(plot_vln_exhaustion_stage_list, function(p) {
            p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
            # stage_palette is now guaranteed to exist in this scope
            p <- p + scale_fill_manual(values = stage_palette, name = "Tumor Stage")
            p # Return the modified plot
        })
        # Combine plots using patchwork
        combined_vln_stage_plot <- wrap_plots(plots = plot_vln_exhaustion_stage, ncol=min(3, length(exhaustion_markers)))
        # Save taller plot
        ggsave(file.path(output_dir, "1_tcell_exhaustion_violin_by_stage.pdf"), combined_vln_stage_plot, width = max(10, length(exhaustion_markers)*3), height = 9)
        cat("Saved: 1_tcell_exhaustion_violin_by_stage.pdf\n")

        # Plot 2: Dot Plot by Stage
        # Add stage to the simplified identity for grouping
        t_cell_subset$celltype_stage <- paste(t_cell_subset$simplified_cell_type, t_cell_subset[[stage_column, drop=TRUE]], sep = "_")
        # Get unique sorted levels for plotting order
        celltype_stage_levels <- sort(unique(t_cell_subset$celltype_stage))
        t_cell_subset$celltype_stage <- factor(t_cell_subset$celltype_stage, levels = celltype_stage_levels)
        
        plot_dot_exhaustion_stage <- DotPlot(t_cell_subset, 
                                             features = exhaustion_markers, 
                                             group.by = "celltype_stage") + # Group by combined identity
                                      RotatedAxis() + 
                                      labs(title = "T Cell Exhaustion Marker Expression by Stage (Simplified Names)",
                                           x = "Cell Type and Stage", y = "Gene") +
                                      theme(axis.text.x = element_text(angle = 45, hjust=1), 
                                            axis.title.x = element_blank()) # Remove x-axis title for space
        ggsave(file.path(output_dir, "2_tcell_exhaustion_dotplot_by_stage.pdf"), plot_dot_exhaustion_stage, width = max(10, length(celltype_stage_levels)*0.5), height = 6) # Adjust width based on number of groups
        cat("Saved: 2_tcell_exhaustion_dotplot_by_stage.pdf\n")

        # --- Add Exhaustion Score Calculation and Plot ---
        cat("Calculating T cell exhaustion module score...\n")
        if (length(exhaustion_markers) > 0) {
            t_cell_subset <- AddModuleScore(t_cell_subset, 
                                            features = list(Exhaustion = exhaustion_markers), 
                                            name = "Exhaustion_Score", 
                                            assay = "RNA", 
                                            seed = 123)
            
            # Plot 3: Violin Plot for Exhaustion Module Score by Stage using ggplot + ggpubr
            score_col_tcell <- "Exhaustion_Score1"
            if (score_col_tcell %in% colnames(t_cell_subset@meta.data)) {
                
                # Define comparisons against "Normal" stage
                if ("Normal" %in% stage_levels && length(stage_levels) >= 2) {
                    other_stages <- stage_levels[stage_levels != "Normal"]
                    my_comparisons <- lapply(other_stages, function(stage) c("Normal", stage))
                } else {
                    my_comparisons <- list()
                    cat("Warning: Cannot perform comparisons against 'Normal' stage. 'Normal' not found or less than 2 stages defined.\n")
                }

                plot_vln_tcell_score_stage <- ggplot(t_cell_subset@meta.data, 
                                                       aes(x = !!sym(stage_column), y = !!sym(score_col_tcell), fill = !!sym(stage_column))) + 
                                           geom_violin(trim = FALSE, scale = "width", draw_quantiles = c(0.5)) + 
                                           facet_wrap(~simplified_cell_type, scales = "free_y") + # Facet by CD4/CD8
                                           scale_fill_manual(values = stage_palette, name = "Tumor Stage") +
                                           labs(title="T Cell Exhaustion Score by Stage", 
                                                x = "Tumor Stage", 
                                                y = "Exhaustion Score") +
                                           theme_bw(base_size = 10) + 
                                           theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                                                 legend.position="right",
                                                 strip.background = element_rect(fill="grey90"))
                
                # Add statistical comparisons if there are stages to compare
                if (length(my_comparisons) > 0) {
                     plot_vln_tcell_score_stage <- plot_vln_tcell_score_stage + 
                        stat_compare_means(comparisons = my_comparisons, 
                                           method = "wilcox.test", 
                                           label = "p.signif", # Show stars: **** <= 0.0001, *** <= 0.001, ** <= 0.01, * <= 0.05
                                           hide.ns = TRUE,       # Hide non-significant comparisons
                                           size = 3)
                } else {
                    cat("Skipping statistical comparisons: Not enough stages defined.\n")
                }
                
                # Increase plot height slightly for significance bars
                ggsave(file.path(output_dir, "3_tcell_exhaustion_score_violin_by_stage.pdf"), plot_vln_tcell_score_stage, width = 8, height = 6)
                cat("Saved: 3_tcell_exhaustion_score_violin_by_stage.pdf\n")
            } else {
                 cat("Skipping T cell exhaustion score violin plot: Score column ('", score_col_tcell, "') not found.\n")
            }
        } else {
            cat("Skipping T cell exhaustion score calculation: No exhaustion markers found.\n")
        }
        # --- End Exhaustion Score --- 

    } else {
        cat("Skipping stage-grouped T cell exhaustion plots because stage column ('", stage_column, "') was not found.\n")
    }

} else {
    cat("Skipping T cell exhaustion analysis as no relevant markers were found in the object.\n")
}

cat("--- Finished T cell exhaustion analysis ---\n")

# --- 2. NK Cell Marker Analysis ---
cat("\n--- Starting NK cell marker analysis ---\n")

# Note: Analyzing 'CD8+ NKT-like cells' and 'ISG expressing immune cells' for NK markers as no specific NK cell type was identified.
if (length(nk_markers) > 0 && length(nk_types) > 0) {
    # Subset Seurat object for NK-like and ISG cells
    nk_subset <- subset(hgsoc_seurat, subset = !!sym(metadata_column) %in% nk_types)

    if (nrow(nk_subset) > 0) {
        # Only generate plot if stage column is available
        if (!is.null(stage_column)) {
            cat("Generating NK marker DotPlot (for NK-like and ISG cells) grouped by stage...\n")
            # Add stage to the identity for grouping
            nk_subset$celltype_stage <- paste(nk_subset[[metadata_column, drop=TRUE]], nk_subset[[stage_column, drop=TRUE]], sep = "_")
            # Get unique sorted levels for plotting order
            nk_celltype_stage_levels <- sort(unique(nk_subset$celltype_stage))
            nk_subset$celltype_stage <- factor(nk_subset$celltype_stage, levels = nk_celltype_stage_levels)

            plot_dot_nk_stage <- DotPlot(nk_subset, 
                                         features = nk_markers, 
                                         group.by = "celltype_stage") + # Group by combined identity
                                  RotatedAxis() + 
                                  labs(title = "NK Marker Expression in NK-like/ISG Cells by Stage",
                                       x = "Cell Type and Stage", y = "Gene") +
                                  theme(axis.text.x = element_text(angle = 45, hjust=1), 
                                        axis.title.x = element_blank()) 
            # Adjust width based on number of groups
            ggsave(file.path(output_dir, "4_nk_marker_dotplot_by_stage.pdf"), plot_dot_nk_stage, width = max(8, length(nk_celltype_stage_levels)*0.6), height = 5)
            cat("Saved: 4_nk_marker_dotplot_by_stage.pdf\n")
        } else {
             cat("Skipping stage-grouped NK marker plot because stage column ('", stage_column, "') was not found. Generating basic plot.\n")
             # Basic DotPlot if no stage info
             plot_dot_nk <- DotPlot(nk_subset, 
                                    features = nk_markers, 
                                    group.by = metadata_column) + 
                              RotatedAxis() + 
                              labs(title = "NK Marker Expression in NK-like/ISG Cells") +
                              theme(axis.text.x = element_text(angle = 45, hjust=1))
             ggsave(file.path(output_dir, "4_nk_marker_dotplot.pdf"), plot_dot_nk, width = 7, height = 5)
             cat("Saved: 4_nk_marker_dotplot.pdf\n")
        }
    } else {
        cat("Skipping NK marker analysis as no cells were found for the specified nk_types: ", paste(nk_types, collapse=", "), "\n")
    }

} else {
    cat("Skipping NK marker analysis as no relevant markers were found or nk_types are not defined.\n")
}

cat("--- Finished NK cell marker analysis ---\n")

# --- 3. M1/M2 Macrophage Score ---
cat("\n--- Starting M1/M2 Macrophage scoring ---\n")

# Check if genes and cell types are valid
if (length(unlist(m1_genes)) > 0 && length(unlist(m2_genes)) > 0 && length(macrophage_types) > 0) {
    
    # Subset for Macrophages
    mac_subset <- subset(hgsoc_seurat, subset = !!sym(metadata_column) %in% macrophage_types)

    if (nrow(mac_subset) > 0) {
        cat("Calculating M1/M2 module scores...\n")
        # Calculate M1 score
        if (length(m1_genes[[1]]) > 0) {
             mac_subset <- AddModuleScore(mac_subset, features = m1_genes, name = "M1_Score", assay = "RNA", seed = 123)
        } else {
            cat("Skipping M1 score calculation: No M1 genes found in the dataset.\n")
            mac_subset$M1_Score1 <- NA # Add placeholder column if score cannot be calculated
        }
        # Calculate M2 score
        if (length(m2_genes[[1]]) > 0) {
            mac_subset <- AddModuleScore(mac_subset, features = m2_genes, name = "M2_Score", assay = "RNA", seed = 123)
        } else {
             cat("Skipping M2 score calculation: No M2 genes found in the dataset.\n")
             mac_subset$M2_Score1 <- NA # Add placeholder column if score cannot be calculated
        }
        
        # Combine individual marker genes for plotting
        individual_m1_m2_markers <- unique(c(m1_genes[[1]], m2_genes[[1]]))
        
        # Generate plots only if stage column is available
        if (!is.null(stage_column)) {
            cat("Generating Macrophage marker and score plots grouped by stage...\n")
            
            # Plot 4: DotPlot for individual M1/M2 markers by stage
            if(length(individual_m1_m2_markers) > 0) {
                mac_subset$celltype_stage <- paste(mac_subset[[metadata_column, drop=TRUE]], mac_subset[[stage_column, drop=TRUE]], sep = "_")
                mac_celltype_stage_levels <- sort(unique(mac_subset$celltype_stage))
                mac_subset$celltype_stage <- factor(mac_subset$celltype_stage, levels = mac_celltype_stage_levels)

                plot_dot_mac_markers_stage <- DotPlot(mac_subset, 
                                                    features = individual_m1_m2_markers, 
                                                    group.by = "celltype_stage") + 
                                          RotatedAxis() + 
                                          labs(title = "M1/M2 Marker Expression in Macrophages by Stage",
                                               x = "Cell Type and Stage", y = "Gene") +
                                          theme(axis.text.x = element_text(angle = 45, hjust=1), 
                                                axis.title.x = element_blank()) 
                ggsave(file.path(output_dir, "5_mac_marker_dotplot_by_stage.pdf"), plot_dot_mac_markers_stage, width = max(8, length(mac_celltype_stage_levels)*0.6), height = max(5, length(individual_m1_m2_markers)*0.3))
                cat("Saved: 5_mac_marker_dotplot_by_stage.pdf\n")
            } else {
                 cat("Skipping Macrophage individual marker DotPlot: No M1/M2 genes found.\n")
            }

            # Plot 5a: Violin plots for individual M1 markers by stage
            if (length(m1_genes[[1]]) > 0) {
                plot_vln_m1_marker_list <- VlnPlot(mac_subset, 
                                                   features = m1_genes[[1]], 
                                                   group.by = stage_column,
                                                   pt.size = 0, combine = FALSE, log = FALSE, raster = FALSE)
                plot_vln_m1_markers <- lapply(plot_vln_m1_marker_list, function(p) p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") + xlab("Tumor Stage"))
                combined_vln_m1_markers <- wrap_plots(plots = plot_vln_m1_markers, ncol=min(3, length(m1_genes[[1]])))
                ggsave(file.path(output_dir, "6a_mac_M1_marker_violin_by_stage.pdf"), combined_vln_m1_markers, width = max(8, length(m1_genes[[1]])*2.5), height = 5)
                cat("Saved: 6a_mac_M1_marker_violin_by_stage.pdf\n")
            } else {
                cat("Skipping M1 individual marker violin plots: No M1 genes found.\n")
            }
            
            # Plot 5b: Violin plots for individual M2 markers by stage
             if (length(m2_genes[[1]]) > 0) {
                plot_vln_m2_marker_list <- VlnPlot(mac_subset, 
                                                   features = m2_genes[[1]], 
                                                   group.by = stage_column,
                                                   pt.size = 0, combine = FALSE, log = FALSE, raster = FALSE)
                plot_vln_m2_markers <- lapply(plot_vln_m2_marker_list, function(p) p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") + xlab("Tumor Stage"))
                combined_vln_m2_markers <- wrap_plots(plots = plot_vln_m2_markers, ncol=min(3, length(m2_genes[[1]])))
                ggsave(file.path(output_dir, "6b_mac_M2_marker_violin_by_stage.pdf"), combined_vln_m2_markers, width = max(8, length(m2_genes[[1]])*2.5), height = 5)
                cat("Saved: 6b_mac_M2_marker_violin_by_stage.pdf\n")
            } else {
                cat("Skipping M2 individual marker violin plots: No M2 genes found.\n")
            }
            
            # Plot 7: Violin plot for M1/M2 Module Scores by stage
            score_cols <- c("M1_Score1", "M2_Score1")[c("M1_Score1", "M2_Score1") %in% colnames(mac_subset@meta.data)] # Check which scores were calculated
            if(length(score_cols) > 0) {
                
                # --- Add Comparisons vs Normal --- 
                # Define comparisons against "Normal" stage
                # Also define palette here for this specific plot
                if ("Normal" %in% stage_levels && length(stage_levels) >= 2) {
                    other_stages <- stage_levels[stage_levels != "Normal"]
                    my_comparisons_mac <- lapply(other_stages, function(stage) c("Normal", stage))
                    
                    # Define YlOrRd palette for macrophage plot stages
                    n_stages <- length(stage_levels)
                    stage_palette_mac <- RColorBrewer::brewer.pal(max(3, n_stages), "YlOrRd")[1:n_stages]
                    names(stage_palette_mac) <- stage_levels
                    
                } else {
                    my_comparisons_mac <- list()
                    stage_palette_mac <- NULL # No palette if stages invalid for comparison
                    cat("Warning: Cannot perform comparisons against 'Normal' stage for Macrophages. 'Normal' not found or less than 2 stages defined.\n")
                }
                # --- End Comparisons --- 

                # Need to pivot longer for ggplot and stats
                mac_scores_long <- mac_subset@meta.data %>% 
                    select(!!sym(stage_column), all_of(score_cols)) %>% 
                    pivot_longer(cols = all_of(score_cols), names_to = "Score_Type", values_to = "Score_Value") %>%
                    mutate(Score_Type = gsub("1$", "", Score_Type)) # Clean up names (e.g., M1_Score1 -> M1_Score)
                
                plot_vln_scores_stage <- ggplot(mac_scores_long, aes(x = !!sym(stage_column), y = Score_Value, fill = !!sym(stage_column))) + 
                                           geom_violin(trim = FALSE, scale = "width", draw_quantiles = c(0.5)) + 
                                           facet_wrap(~Score_Type, scales = "free_y") + # Facet by M1/M2 score
                                           # Use the specific macrophage stage palette
                                           {if (!is.null(stage_palette_mac)) scale_fill_manual(values = stage_palette_mac, name = "Tumor Stage") else NULL} +
                                           labs(title="Macrophage M1 vs M2 Scores by Stage", 
                                                x = "Tumor Stage", 
                                                y = "Module Score") +
                                           theme_bw(base_size = 10) + 
                                           theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                                                 legend.position="right",
                                                 strip.background = element_rect(fill="grey90"))
                
                # Add statistical comparisons if defined
                if (length(my_comparisons_mac) > 0) {
                    plot_vln_scores_stage <- plot_vln_scores_stage + 
                        stat_compare_means(comparisons = my_comparisons_mac, 
                                           method = "wilcox.test", 
                                           label = "p.signif", 
                                           hide.ns = TRUE, 
                                           size = 3)
                } 
                
                # Adjust plot height for facets and stats
                ggsave(file.path(output_dir, "7_mac_module_score_violin_by_stage.pdf"), plot_vln_scores_stage, width = 8, height = 6)
                cat("Saved: 7_mac_module_score_violin_by_stage.pdf\n")
            } else {
                cat("Skipping module score violin plot: Neither M1 nor M2 scores could be calculated.\n")
            }

        } else {
             cat("Skipping stage-grouped Macrophage plots because stage column ('", stage_column, "') was not found. Consider adding basic plots if needed.\n")
             # Add non-stage-grouped plots here if desired
        }
    } else {
        cat("Skipping M1/M2 analysis as no cells were found for the specified macrophage_types: ", paste(macrophage_types, collapse=", "), "\n")
    }
} else {
    cat("Skipping M1/M2 analysis as no relevant M1/M2 markers were found or macrophage_types are not defined.\n")
}

cat("--- Finished M1/M2 Macrophage scoring ---\n")

cat("\nFunctional analysis script finished.\n") 