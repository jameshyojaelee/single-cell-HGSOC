# =====================================================
# Gene Set Enrichment Analysis (GSEA) of Non-Immune DEGs
# =====================================================

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db) # Human organism database
library(dplyr)
library(readr)
library(stringr)
library(tibble)

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/diff_expr/non_immune")
output_dir <- file.path(base_dir, "output/enrichment/GSEA")
plot_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results")

# Create output directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Parameters
gsea_ontology <- "BP" # Using GO Biological Process for GSEA
num_gsea_plots <- 5 # Number of top pathways to show in gseaplot2

# --- Helper Function for GSEA Analysis ---

run_gsea_analysis <- function(deg_file_path, comparison_name) {
  cat(paste("\nProcessing GSEA for:", comparison_name, "\n"))
  
  # Read DEG results
  deg_results <- read_csv(deg_file_path, show_col_types = FALSE)
  
  # Prepare ranked gene list
  # 1. Filter out genes with NA logFC values
  # 2. Handle potential duplicate gene symbols (take the one with max abs(logFC))
  # 3. Convert to Entrez IDs
  # 4. Create named vector sorted by logFC
  
  ranked_genes <- deg_results %>%
    filter(!is.na(avg_log2FC)) %>%
    group_by(gene) %>%
    filter(abs(avg_log2FC) == max(abs(avg_log2FC))) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE) # Keep only one entry per gene symbol
  
  if (nrow(ranked_genes) == 0) {
      cat("  No genes with valid logFC found. Skipping GSEA.\n")
      return(NULL)
  }

  # Convert gene symbols to Entrez IDs
  entrez_map <- tryCatch({
    bitr(ranked_genes$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  }, error = function(e) {
      cat("  Error during Gene ID conversion: ", e$message, "\n")
      return(NULL)
  })
  
  if (is.null(entrez_map) || nrow(entrez_map) == 0) {
    cat("  No Entrez IDs could be mapped. Skipping GSEA.\n")
    return(NULL)
  }
  
  # Merge Entrez IDs back and prepare the final ranked list
  ranked_list_df <- ranked_genes %>%
    inner_join(entrez_map, by = c("gene" = "SYMBOL")) %>%
    select(ENTREZID, avg_log2FC) %>%
    # Handle potential duplicate Entrez IDs if multiple symbols map to one ID
    group_by(ENTREZID) %>%
    summarise(avg_log2FC = mean(avg_log2FC)) %>%
    ungroup()

  # Create the final named vector, sorted by logFC
  gene_list <- ranked_list_df$avg_log2FC
  names(gene_list) <- ranked_list_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat(paste("  Prepared ranked list with", length(gene_list), "genes.\n"))
  
  # Run GSEA using gseGO
  cat(paste("  Running gseGO for Ontology:", gsea_ontology, "...\n"))
  
  tryCatch({
    gsea_result <- gseGO(
      geneList     = gene_list,
      OrgDb        = org.Hs.eg.db,
      ont          = gsea_ontology,
      minGSSize    = 10, # Minimum size of gene set for analysis
      maxGSSize    = 500, # Maximum size of gene set for analysis
      pvalueCutoff = 0.05,
      verbose      = FALSE, # Set to TRUE for more detailed output
      pAdjustMethod= "BH"
    )
    
    if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
      cat(paste("  No significant GSEA enrichment found for", gsea_ontology, ".\n"))
      return(NULL)
    }
    
    # Save GSEA results table
    gsea_df <- as.data.frame(gsea_result@result)
    write.csv(gsea_df, 
              file.path(results_dir, paste0("GSEA_", comparison_name, "_", gsea_ontology, "_results.csv")),
              row.names = FALSE)
    cat(paste("  Saved GSEA table for", gsea_ontology, ".\n"))
    
    # Generate and save plots
    # GSEA Plot for top N pathways
    plot_indices <- head(order(abs(gsea_result$NES), decreasing = TRUE), num_gsea_plots)
    p_gsea <- gseaplot2(gsea_result, geneSetID = plot_indices, title = paste("GSEA Plot:", comparison_name))
    ggsave(file.path(plot_dir, paste0("GSEA_", comparison_name, "_", gsea_ontology, "_gseaplot.pdf")),
            p_gsea, width = 10, height = 8)
            
    # Ridge Plot
    p_ridge <- ridgeplot(gsea_result, showCategory = num_gsea_plots) + 
                 ggtitle(paste("GSEA Ridge Plot:", comparison_name))
    ggsave(file.path(plot_dir, paste0("GSEA_", comparison_name, "_", gsea_ontology, "_ridgeplot.pdf")),
            p_ridge, width = 10, height = 8)

    cat(paste("  Saved GSEA plots for", gsea_ontology, ".\n"))
      
  }, error = function(e) {
    cat(paste("  Error during GSEA analysis for", gsea_ontology, ":", e$message, "\n"))
  })
}

# --- Main Execution --- 

# Find all DEG result files
deg_files <- list.files(path = input_dir, pattern = "^DEG_.*_vs_Normal_wilcox\\.csv$", full.names = TRUE)

if (length(deg_files) == 0) {
  stop(paste("No DEG result files found in:", input_dir))
}

cat(paste("Found", length(deg_files), "DEG files to process for GSEA.\n"))

# Process each file
for (file_path in deg_files) {
  # Extract comparison name from filename (e.g., IC2_vs_Normal)
  comparison_name <- basename(file_path) %>%
    str_replace("^DEG_", "") %>%
    str_replace("_wilcox\\.csv$", "") # Corrected regex escape
  
  run_gsea_analysis(file_path, comparison_name)
}

cat("\nGSEA analysis complete.\n")
cat(paste("Results saved in:", output_dir, "\n")) 