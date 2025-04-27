# ==========================================================================
# Transcription Factor (TF) Enrichment Analysis using RcisTarget
# ==========================================================================

# Load required libraries
library(RcisTarget)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)

# --- Configuration ---

# Path to the directory containing the downloaded RcisTarget feather databases
db_path <- "/gpfs/commons/home/jameslee/HGSOC/databases/RcisTarget/hg38" # <<< UPDATE THIS PATH

# Check if db_path exists
if (!dir.exists(db_path)) {
  stop(paste("RcisTarget database directory not found:", db_path, 
             "\nPlease download the databases and update the db_path variable."))
}

# Specify database files (adjust filenames if necessary based on your downloads)
motif_rankings_file <- file.path(db_path, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather") 
motif_annotations_file <- file.path(db_path, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

# Check if database files exist
if (!file.exists(motif_rankings_file)) stop(paste("Motif rankings file not found:", motif_rankings_file))
if (!file.exists(motif_annotations_file)) stop(paste("Motif annotations file not found:", motif_annotations_file))

# Other paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
input_dir <- file.path(base_dir, "output/diff_expr/non_immune")
output_dir <- file.path(base_dir, "output/enrichment/TF")
plot_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results")

# Create output directories if they don't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Parameters
p_adj_threshold <- 0.05
logfc_threshold <- 0.25 # Consistent with DEG script? Check if needed.
num_tfs_plot <- 15 # Number of top TFs to plot

# --- RcisTarget Initialization ---

cat("Loading RcisTarget databases...\n")
motifRankings <- importRankings(motif_rankings_file)
motifAnnot <- read.table(motif_annotations_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)
cat("Databases loaded.\n")

# --- Helper Function for TF Analysis ---

run_tf_enrichment <- function(deg_file_path, comparison_name) {
  cat(paste("\nProcessing TF enrichment for:", comparison_name, "\n"))
  
  deg_results <- read_csv(deg_file_path, show_col_types = FALSE)
  
  # Separate up- and down-regulated genes
  genes_up <- deg_results %>%
    filter(p_val_adj < p_adj_threshold & avg_log2FC > logfc_threshold) %>%
    pull(gene)
  genes_down <- deg_results %>%
    filter(p_val_adj < p_adj_threshold & avg_log2FC < -logfc_threshold) %>%
    pull(gene)
  
  process_gene_set <- function(gene_set, direction) {
    if (length(gene_set) < 5) { # Need a minimum number of genes
      cat(paste("  Skipping", direction, "regulated: Too few significant genes (", length(gene_set), ").\n"))
      return(NULL)
    }
    cat(paste("  Analyzing", length(gene_set), direction, "regulated genes...\n"))
    
    tryCatch({
      # Calculate enrichment
      motifs_AUC <- calcAUC(gene_set, motifRankings, nCores=4) # Adjust nCores
      motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                               motifAnnot=motifAnnot, 
                                               motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                                               motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity", 
                                                                   "inferredBy_MotifSimilarity_n_Orthology"))
      
      # Filter for significant motifs/TFs (High confidence annotation & significant AUC)
      enriched_tfs <- motifEnrichmentTable %>%
        filter(AUC > 0.01) # Filter based on AUC threshold (adjust as needed)
      
      if (nrow(enriched_tfs) == 0) {
          cat(paste("  No significantly enriched TFs found for", direction, "regulated genes.\n"))
          return(NULL)
      }

      # Keep relevant columns and add TF name
      enriched_tfs <- enriched_tfs %>%
        select(motif, AUC, NES, highConfAnnot, TF_highConf) %>%
        rename(TF = TF_highConf) %>% # Use high confidence TF annotations
        filter(!is.na(TF) & TF != "") %>%
        arrange(desc(NES))

      # Save results table
      results_file <- file.path(results_dir, paste0("TF_", comparison_name, "_", direction, "_enrichment.csv"))
      write.csv(enriched_tfs, results_file, row.names = FALSE)
      cat(paste("  Saved TF enrichment table for", direction, "regulated genes.\n"))
      
      # Generate Bar Plot for top TFs
      top_tfs_plot <- head(enriched_tfs, num_tfs_plot)
      
      p_bar <- ggplot(top_tfs_plot, aes(x = reorder(TF, NES), y = NES, fill = NES)) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient(low = "lightblue", high = "darkblue", name = "NES") +
        labs(
          title = paste("Top Enriched TFs:", comparison_name, "(", direction, "reg.)"),
          x = "Transcription Factor",
          y = "Normalized Enrichment Score (NES)"
        ) +
        theme_minimal(base_size = 10)
        
      plot_file <- file.path(plot_dir, paste0("TF_", comparison_name, "_", direction, "_barplot.pdf"))
      ggsave(plot_file, p_bar, width = 8, height = 7)
      cat(paste("  Saved TF enrichment plot for", direction, "regulated genes.\n"))
      
      return(enriched_tfs)
      
    }, error = function(e) {
      cat(paste("  Error during RcisTarget analysis for", direction, "regulated genes:", e$message, "\n"))
      return(NULL)
    })
  }
  
  # Process both sets
  results_up <- process_gene_set(genes_up, "up")
  results_down <- process_gene_set(genes_down, "down")
  
  # Optional: Combine results or further analysis here if needed
}

# --- Main Execution --- 

# Find all DEG result files
deg_files <- list.files(path = input_dir, pattern = "^DEG_.*_vs_Normal_wilcox\\.csv$", full.names = TRUE)

if (length(deg_files) == 0) {
  stop(paste("No DEG result files found in:", input_dir))
}

cat(paste("Found", length(deg_files), "DEG files to process for TF analysis.\n"))

# Process each file
for (file_path in deg_files) {
  comparison_name <- basename(file_path) %>%
    str_replace("^DEG_", "") %>%
    str_replace("_wilcox\\.csv$", "") # Corrected regex escape
  
  run_tf_enrichment(file_path, comparison_name)
}

cat("\nTF enrichment analysis complete.\n")
cat(paste("Results saved in:", output_dir, "\n")) 