# ============================================================
# Gene Ontology (GO) Enrichment Analysis of Non-Immune DEGs
# ============================================================

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db) # Human organism database
library(dplyr)
library(readr)
library(stringr)

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/diff_expr/non_immune")
output_dir <- file.path(base_dir, "output/enrichment/GO")
plot_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results")

# Create output directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Parameters
p_adj_threshold <- 0.05
go_ontologies <- c("BP", "MF", "CC") # Biological Process, Molecular Function, Cellular Component
num_terms_plot <- 15 # Number of top terms to show in plots

# --- Helper Function for GO Analysis ---

run_go_analysis <- function(deg_file_path, comparison_name) {
  cat(paste("\nProcessing:", comparison_name, "\n"))
  
  # Read DEG results
  deg_results <- read_csv(deg_file_path, show_col_types = FALSE)
  
  # Filter significant genes
  sig_genes <- deg_results %>%
    filter(p_val_adj < p_adj_threshold) %>%
    pull(gene)
  
  if (length(sig_genes) == 0) {
    cat("  No significant genes found. Skipping GO analysis.\n")
    return(NULL)
  }
  
  cat(paste("  Found", length(sig_genes), "significant genes.\n"))
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- tryCatch({
    bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  }, error = function(e) {
      cat("  Error during Gene ID conversion: ", e$message, "\n")
      return(NULL)
  })
  
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    cat("  No Entrez IDs could be mapped. Skipping GO analysis.\n")
    return(NULL)
  }
  
  cat(paste("  Mapped", nrow(entrez_ids), "genes to Entrez IDs.\n"))
  
  # Run enrichGO for each ontology
  for (ont in go_ontologies) {
    cat(paste("  Running enrichGO for Ontology:", ont, "...\n"))
    
    tryCatch({
      ego <- enrichGO(
        gene          = entrez_ids$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        keyType       = 'ENTREZID',
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05, # Standard cutoff for enrichment p-value
        qvalueCutoff  = 0.2   # Standard cutoff for enrichment q-value
      )
      
      if (is.null(ego) || nrow(ego@result) == 0) {
        cat(paste("  No significant enrichment found for", ont, ".\n"))
        next
      }
      
      # Simplify results (optional, removes redundancy)
      ego_simple <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      
      if (nrow(ego_simple@result) == 0) {
        cat(paste("  No significant enrichment found after simplify() for", ont, ". Using original results.\n"))
        ego_to_plot <- ego
      } else {
           cat(paste("  Simplified enrichment results for", ont, ".\n"))
        ego_to_plot <- ego_simple
      }
       
      # Save simplified results table
      results_df <- as.data.frame(ego_to_plot@result)
      write.csv(results_df, 
                file.path(results_dir, paste0("GO_", comparison_name, "_", ont, "_enrichment.csv")),
                row.names = FALSE)
      cat(paste("  Saved enrichment table for", ont, ".\n"))
      
      # Generate and save plots
      # Dot Plot
      p_dot <- dotplot(ego_to_plot, showCategory=num_terms_plot) + 
                ggtitle(paste("GO", ont, "Dot Plot:", comparison_name))
      ggsave(file.path(plot_dir, paste0("GO_", comparison_name, "_", ont, "_dotplot.pdf")),
             p_dot, width = 10, height = 8)
             
      # Bar Plot
      p_bar <- barplot(ego_to_plot, showCategory=num_terms_plot) +
               ggtitle(paste("GO", ont, "Bar Plot:", comparison_name))
      ggsave(file.path(plot_dir, paste0("GO_", comparison_name, "_", ont, "_barplot.pdf")),
              p_bar, width = 10, height = 8)

      cat(paste("  Saved plots for", ont, ".\n"))
      
    }, error = function(e) {
      cat(paste("  Error during enrichGO for", ont, ":", e$message, "\n"))
    })
  }
}

# --- Main Execution --- 

# Find all DEG result files
deg_files <- list.files(path = input_dir, pattern = "^DEG_.*_vs_Normal_wilcox\\.csv$", full.names = TRUE)

if (length(deg_files) == 0) {
  stop(paste("No DEG result files found in:", input_dir))
}

cat(paste("Found", length(deg_files), "DEG files to process.\n"))

# Process each file
for (file_path in deg_files) {
  # Extract comparison name from filename (e.g., IC2_vs_Normal)
  comparison_name <- basename(file_path) %>%
    str_replace("^DEG_", "") %>%
    str_replace("_wilcox\\\\.csv$", "") # Escaping backslash and dot
  
  run_go_analysis(file_path, comparison_name)
}

cat("\nGO enrichment analysis complete.\n")
cat(paste("Results saved in:", output_dir, "\n")) 