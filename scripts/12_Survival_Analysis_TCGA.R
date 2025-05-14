# ============================================================
# Survival Analysis of DEGs using TCGA-OV Data
# ============================================================

# Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)     # Add tibble for rownames_to_column function
library(TCGAbiolinks)
library(glmnet)     # For LASSO/Ridge regression
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)
library(SummarizedExperiment)
library(fgsea)      # For gene set enrichment

# Set error handling to get more information
options(error = function() {
  traceback(2)
  if(!interactive()) quit(status = 1)
})

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
deg_dir <- file.path(base_dir, "output/diff_expr")
output_dir <- file.path(base_dir, "output/survival_analysis")
plot_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results")

# Create output directories if they don't exist
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# List of cell types to analyze
cell_types <- c("non_immune", "t_cells", "macrophages")

# Parameters
p_value_threshold <- 0.05
q_value_threshold <- 0.1
top_n_genes <- 20     # Number of top genes to visualize

# --- Download and Process TCGA-OV Data ---

cat("Downloading TCGA-OV clinical data...\n")

# Use GDCquery_clinic instead of the BCR Biotab approach
clinical_data <- GDCquery_clinic(project = "TCGA-OV", type = "clinical")

cat("Clinical data downloaded:", nrow(clinical_data), "patients\n")

# First, inspect the column names to understand what we're working with
cat("Available clinical data columns:\n")
print(colnames(clinical_data))

# Save the clinical data for reference
write.csv(head(clinical_data, 10), file.path(output_dir, "clinical_data_sample.csv"))

# Standardize column names - with checks to ensure columns exist
clinical_standardized <- clinical_data

# Identify key columns needed for survival analysis
id_col <- NULL
vital_status_col <- NULL
survival_time_col <- NULL
death_time_col <- NULL

# Find patient ID column (usually 'bcr_patient_barcode' or 'submitter_id')
id_candidates <- c("bcr_patient_barcode", "submitter_id", "case_submitter_id", "patient_id")
for (col in id_candidates) {
  if (col %in% colnames(clinical_data)) {
    id_col <- col
    break
  }
}
if (is.null(id_col)) {
  stop("Could not find patient ID column in clinical data")
}

# Find vital status column
status_candidates <- c("vital_status", "patient.vital_status", "OS_status")
for (col in status_candidates) {
  if (col %in% colnames(clinical_data)) {
    vital_status_col <- col
    break
  }
}
if (is.null(vital_status_col)) {
  stop("Could not find vital status column in clinical data")
}

# Find survival time columns
time_candidates <- c("days_to_last_follow_up", "OS_time", "survival_time", "days_to_last_followup")
for (col in time_candidates) {
  if (col %in% colnames(clinical_data)) {
    survival_time_col <- col
    break
  }
}

death_candidates <- c("days_to_death", "death_time", "days_to_patient_death")
for (col in death_candidates) {
  if (col %in% colnames(clinical_data)) {
    death_time_col <- col
    break
  }
}

# If we couldn't find at least one time column, we can't do survival analysis
if (is.null(survival_time_col) && is.null(death_time_col)) {
  stop("Could not find survival time or death time column in clinical data")
}

# Create standardized column names
clinical_data$patient_id <- clinical_data[[id_col]]
clinical_data$vital_status <- clinical_data[[vital_status_col]]

# Handle survival time - use days_to_death when available, otherwise use follow-up time
if (!is.null(death_time_col)) {
  clinical_data$days_to_death <- clinical_data[[death_time_col]]
}

if (!is.null(survival_time_col)) {
  clinical_data$overall_survival_time <- clinical_data[[survival_time_col]]
} else {
  clinical_data$overall_survival_time <- NA
}

# Prioritize days_to_death for survival time when it's available
if (!is.null(death_time_col)) {
  for (i in 1:nrow(clinical_data)) {
    if (!is.na(clinical_data$days_to_death[i])) {
      clinical_data$overall_survival_time[i] <- clinical_data$days_to_death[i]
    }
  }
}

# Calculate survival status (1 = dead, 0 = alive)
clinical_data$overall_survival_status <- ifelse(
  tolower(as.character(clinical_data$vital_status)) %in% c("dead", "deceased"), 1, 0)

# Remove patients with missing survival data or zero/negative survival time
clinical_data <- clinical_data %>%
  filter(!is.na(overall_survival_time) & overall_survival_time > 0)

cat("Clinical data processed:", nrow(clinical_data), "patients\n")

# Download RNA-Seq expression data
cat("Downloading TCGA-OV RNA-Seq data (this may take some time)...\n")
query_exp <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp)
expression_se <- GDCprepare(query_exp)

# Get count data and convert to TPM for better comparability with single-cell data
cat("Processing expression data...\n")

# Create a variable to track if we need to use backup approach
use_backup <- FALSE

# Get available assays
cat("Available assays in SummarizedExperiment:", "\n")
avail_assays <- assays(expression_se)
print(names(avail_assays))

# Choose an appropriate assay - prefer unstranded, but fall back to others if needed
assay_to_use <- "unstranded"
if (!"unstranded" %in% names(avail_assays)) {
  # Try tpm_unstrand which is already normalized
  if ("tpm_unstrand" %in% names(avail_assays)) {
    cat("Using pre-computed tpm_unstrand assay instead of calculating TPM manually\n")
    expression_data <- as.data.frame(assay(expression_se, "tpm_unstrand"))
    expression_data$gene_id <- rownames(expression_data)
    use_backup <- TRUE
  } else if (length(names(avail_assays)) > 0) {
    # Use the first available assay
    assay_to_use <- names(avail_assays)[1]
    cat("Using", assay_to_use, "assay instead of unstranded\n")
  }
}

if (!use_backup) {
  # Get count data from chosen assay
  counts <- assay(expression_se, assay_to_use)
  cat("Using", assay_to_use, "assay for counts\n")
  
  # Print dimensions for debugging
  cat("Dimensions of counts:", dim(counts)[1], "x", dim(counts)[2], "\n")
  
  # Check if gene_lengths are available
  if (!("gene_type_length" %in% colnames(rowData(expression_se)))) {
    cat("Warning: gene_type_length not found in rowData. Using raw counts instead of TPM.\n")
    expression_data <- as.data.frame(counts)
    expression_data$gene_id <- rownames(expression_data)
    use_backup <- TRUE
  }
}

if (!use_backup) {
  gene_lengths <- rowData(expression_se)$gene_type_length
  cat("Length of gene_lengths:", length(gene_lengths), "\n")
  
  # Check if gene_lengths has any NA or zero values
  if (any(is.na(gene_lengths)) || any(gene_lengths == 0)) {
    cat("Warning: gene_lengths contains NA or zero values. Replacing with minimum non-zero value.\n")
    min_length <- min(gene_lengths[gene_lengths > 0], na.rm = TRUE)
    gene_lengths[is.na(gene_lengths) | gene_lengths == 0] <- min_length
  }
  
  # Calculate TPM - with error handling
  tryCatch({
    # Ensure counts is a matrix
    if (!is.matrix(counts)) {
      cat("Converting counts to matrix...\n")
      counts <- as.matrix(counts)
    }
    
    # Calculate RPK (reads per kilobase)
    rpk <- counts / (gene_lengths / 1000)
    
    # Check if rpk is a matrix with correct dimensions
    if (!is.matrix(rpk) || nrow(rpk) == 0 || ncol(rpk) == 0) {
      stop("RPK calculation resulted in an invalid matrix")
    }
    
    # Calculate scaling factor
    cat("Calculating scaling factors...\n")
    scaling_factor <- colSums(rpk, na.rm = TRUE) / 1e6
    
    # Check for zero scaling factors
    if (any(scaling_factor == 0)) {
      cat("Warning: Some scaling factors are zero. Replacing with minimum non-zero value.\n")
      min_sf <- min(scaling_factor[scaling_factor > 0])
      scaling_factor[scaling_factor == 0] <- min_sf
    }
    
    # Calculate TPM
    tpm_matrix <- t(t(rpk) / scaling_factor)
    
    # Convert to data frame
    expression_data <- as.data.frame(tpm_matrix)
    expression_data$gene_id <- rownames(expression_data)
    
  }, error = function(e) {
    cat("Error in TPM calculation:", e$message, "\n")
    cat("Attempting alternative approach using raw counts...\n")
    
    # Fallback to use raw counts instead of TPM
    expression_data <<- as.data.frame(counts)
    expression_data$gene_id <<- rownames(expression_data)
  })
}

# Double-check that expression_data exists
if (!exists("expression_data")) {
  cat("Emergency fallback: creating expression_data from counts\n")
  counts <- assay(expression_se, assay_to_use)
  expression_data <- as.data.frame(counts)
  expression_data$gene_id <- rownames(expression_data)
}

# Map Ensembl IDs to gene symbols
cat("Mapping Ensembl IDs to gene symbols...\n")
gene_info <- rowData(expression_se)

# Get gene names - check column names first
gene_name_col <- NULL
possible_cols <- c("gene_name", "external_gene_name", "Symbol", "gene_symbol")
for (col in possible_cols) {
  if (col %in% colnames(gene_info)) {
    gene_name_col <- col
    break
  }
}

if (is.null(gene_name_col)) {
  cat("Could not find gene name column. Using Ensembl IDs as gene symbols.\n")
  gene_id_map <- data.frame(
    gene_id = rownames(gene_info),
    gene_symbol = rownames(gene_info),
    stringsAsFactors = FALSE
  )
} else {
  gene_id_map <- data.frame(
    gene_id = rownames(gene_info),
    gene_symbol = gene_info[[gene_name_col]],
    stringsAsFactors = FALSE
  )
}

# Add gene symbols to expression data
expression_data <- merge(expression_data, gene_id_map, by = "gene_id", all.x = FALSE)
expression_data <- expression_data[!is.na(expression_data$gene_symbol), ]

# Process expressions to handle multiple Ensembl IDs per gene
genes <- unique(expression_data$gene_symbol)
samples <- colnames(expression_data)[!colnames(expression_data) %in% c("gene_id", "gene_symbol")]

# Create a matrix for the processed data
processed_data <- matrix(NA, nrow = length(genes), ncol = length(samples))
rownames(processed_data) <- genes
colnames(processed_data) <- samples

# Fill in the matrix with max expression per gene
cat("Aggregating expression values by gene symbol...\n")
for (g in 1:length(genes)) {
  if (g %% 1000 == 0) cat("Processed", g, "of", length(genes), "genes\n")
  gene_rows <- expression_data$gene_symbol == genes[g]
  for (s in 1:length(samples)) {
    sample_col <- samples[s]
    processed_data[g, s] <- max(expression_data[gene_rows, sample_col], na.rm = TRUE)
  }
}

# Convert processed matrix to data frame
expression_data_final <- as.data.frame(processed_data)

cat("Expression data processed:", nrow(expression_data_final), "genes,", ncol(expression_data_final), "samples\n")

# Extract TCGA patient IDs from column names (TCGA-XX-XXXX)
sample_to_patient <- data.frame(
  sample_id = colnames(expression_data_final),
  patient_id = str_extract(colnames(expression_data_final), "TCGA-[A-Z0-9]+-[A-Z0-9]+"),
  stringsAsFactors = FALSE
)

# --- Univariate Cox Regression Analysis ---

run_univariate_cox <- function(cell_type) {
  cat("\n=== Running Univariate Cox Regression for", cell_type, "===\n")
  
  # Load DEG results
  deg_file_pattern <- paste0("DEG_.*_vs_Normal_wilcox\\.csv$")
  deg_files <- list.files(path = file.path(deg_dir, cell_type), 
                         pattern = deg_file_pattern, 
                         full.names = TRUE)
  
  if (length(deg_files) == 0) {
    cat("No DEG files found for", cell_type, "\n")
    return(NULL)
  }
  
  # Initialize list to store results
  all_results <- list()
  
  # Process each comparison
  for (deg_file in deg_files) {
    comparison_name <- basename(deg_file) %>%
      str_replace("^DEG_", "") %>%
      str_replace("_wilcox\\.csv$", "")
    
    cat("Processing comparison:", comparison_name, "\n")
    
    # Read DEG results
    deg_results <- read_csv(deg_file, show_col_types = FALSE)
    
    # Extract significant DEGs (adjust p-value threshold as needed)
    significant_degs <- deg_results %>%
      filter(p_val_adj < q_value_threshold) %>%
      arrange(p_val_adj)
    
    if (nrow(significant_degs) == 0) {
      cat("No significant DEGs found for", comparison_name, "\n")
      next
    }
    
    cat("Found", nrow(significant_degs), "significant DEGs\n")
    
    # Extract genes for TCGA data
    deg_genes <- significant_degs$gene
    
    # Subset TCGA expression data for these genes
    # Note: Some genes might not be found in the TCGA data (symbol mismatches, etc.)
    genes_in_expr <- rownames(expression_data_final)[rownames(expression_data_final) %in% deg_genes]
    
    if (length(genes_in_expr) == 0) {
      cat("None of the DEGs found in TCGA expression data for", comparison_name, "\n")
      next
    }
    
    cat("Found", length(genes_in_expr), "DEGs in TCGA expression data\n")
    
    # Prepare data for survival analysis
    gene_expr_subset <- expression_data_final[genes_in_expr, , drop = FALSE]
    
    # Transpose for patient-level analysis (patients as rows)
    gene_expr_subset <- t(gene_expr_subset)
    gene_expr_subset <- as.data.frame(gene_expr_subset) %>%
      rownames_to_column("sample_id")
    
    # Add patient IDs
    gene_expr_subset <- gene_expr_subset %>%
      left_join(sample_to_patient, by = "sample_id") %>%
      select(-sample_id)
    
    # Merge with clinical data
    merged_data <- gene_expr_subset %>%
      inner_join(clinical_data, by = "patient_id")
    
    cat("Merged data:", nrow(merged_data), "patients with both expression and clinical data\n")
    
    # Run univariate Cox for each gene
    cox_results <- data.frame()
    
    for (gene in genes_in_expr) {
      tryCatch({
        # Create Cox formula
        formula <- as.formula(paste("Surv(overall_survival_time, overall_survival_status) ~", gene))
        
        # Fit Cox model
        cox_model <- coxph(formula, data = merged_data)
        
        # Extract results
        cox_summary <- summary(cox_model)
        
        # Store results
        gene_result <- data.frame(
          gene = gene,
          hazard_ratio = exp(cox_summary$coefficients[1, "coef"]),
          ci_lower = exp(cox_summary$coefficients[1, "coef"] - 1.96 * cox_summary$coefficients[1, "se(coef)"]),
          ci_upper = exp(cox_summary$coefficients[1, "coef"] + 1.96 * cox_summary$coefficients[1, "se(coef)"]),
          p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
          significance = ifelse(cox_summary$coefficients[1, "Pr(>|z|)"] < p_value_threshold, 
                                "Significant", "Not Significant"),
          log_rank = cox_summary$logtest[1],
          comparison = comparison_name
        )
        
        cox_results <- rbind(cox_results, gene_result)
        
      }, error = function(e) {
        cat("Error in Cox regression for gene", gene, ":", e$message, "\n")
      })
    }
    
    # Adjust p-values
    if (nrow(cox_results) > 0) {
      cox_results$fdr <- p.adjust(cox_results$p_value, method = "BH")
      
      # Sort by p-value
      cox_results <- cox_results %>%
        arrange(p_value)
      
      # Save results
      results_file <- file.path(results_dir, paste0("univariate_cox_", cell_type, "_", comparison_name, ".csv"))
      write.csv(cox_results, results_file, row.names = FALSE)
      cat("Saved univariate Cox results to", results_file, "\n")
      
      # Store for combined analysis
      all_results[[comparison_name]] <- cox_results
      
      # Visualize top significant genes
      if (any(cox_results$p_value < p_value_threshold)) {
        cat("Creating forest plot for significant genes...\n")
        
        # Select top genes
        top_genes <- cox_results %>%
          filter(p_value < p_value_threshold) %>%
          head(top_n_genes)
        
        # Forest plot
        forest_plot <- ggplot(top_genes, aes(x = reorder(gene, hazard_ratio), y = hazard_ratio)) +
          geom_point(aes(color = significance), size = 3) +
          geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = significance), width = 0.2) +
          geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
          scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "grey")) +
          coord_flip() +
          labs(
            title = paste("Top Survival-Associated DEGs -", comparison_name),
            x = "Gene",
            y = "Hazard Ratio (95% CI)",
            subtitle = paste("Cell type:", cell_type)
          ) +
          theme_bw()
        
        # Save plot
        forest_plot_file <- file.path(plot_dir, paste0("forest_plot_", cell_type, "_", comparison_name, ".pdf"))
        ggsave(forest_plot_file, forest_plot, width = 10, height = 8)
        cat("Saved forest plot to", forest_plot_file, "\n")
        
        # For top genes, create K-M plots (top 5)
        top_5_genes <- top_genes$gene[1:min(5, nrow(top_genes))]
        
        for (gene in top_5_genes) {
          cat("Creating Kaplan-Meier plot for gene", gene, "...\n")
          
          # Create high/low expression groups based on median
          merged_data$expression_group <- ifelse(
            merged_data[[gene]] > median(merged_data[[gene]]), 
            "High", "Low"
          )
          
          # Fit survival curve
          km_fit <- survfit(Surv(overall_survival_time, overall_survival_status) ~ expression_group, 
                           data = merged_data)
          
          # Plot K-M curve
          km_plot <- ggsurvplot(
            km_fit,
            data = merged_data,
            risk.table = TRUE,
            pval = TRUE,
            conf.int = TRUE,
            title = paste("Kaplan-Meier Curve for", gene),
            ggtheme = theme_bw(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE,
            palette = c("blue", "red")
          )
          
          # Save plot
          km_plot_file <- file.path(plot_dir, paste0("km_plot_", cell_type, "_", comparison_name, "_", gene, ".pdf"))
          pdf(km_plot_file, width = 10, height = 8)
          print(km_plot)
          dev.off()
          cat("Saved K-M plot to", km_plot_file, "\n")
        }
      }
    } else {
      cat("No significant Cox regression results for", comparison_name, "\n")
    }
  }
  
  return(all_results)
}

# --- Multivariate and LASSO Cox Regression ---

run_multivariate_cox <- function(cell_type, univariate_results) {
  cat("\n=== Running Multivariate/LASSO Cox Regression for", cell_type, "===\n")
  
  if (is.null(univariate_results) || length(univariate_results) == 0) {
    cat("No univariate results available. Skipping multivariate analysis.\n")
    return(NULL)
  }
  
  # Process each comparison
  for (comparison_name in names(univariate_results)) {
    cat("Processing comparison:", comparison_name, "\n")
    
    # Get univariate results
    uni_results <- univariate_results[[comparison_name]]
    
    # Select significant genes
    significant_genes <- uni_results %>%
      filter(p_value < p_value_threshold) %>%
      pull(gene)
    
    if (length(significant_genes) == 0) {
      cat("No significant genes for multivariate analysis in", comparison_name, "\n")
      next
    }
    
    if (length(significant_genes) < 2) {
      cat("Only one significant gene. Skipping multivariate analysis for", comparison_name, "\n")
      next
    }
    
    cat("Found", length(significant_genes), "significant genes for multivariate analysis\n")
    
    # Subset TCGA expression data
    gene_expr_subset <- expression_data_final[significant_genes, , drop = FALSE]
    gene_expr_subset <- t(gene_expr_subset)
    gene_expr_subset <- as.data.frame(gene_expr_subset) %>%
      rownames_to_column("sample_id")
    
    # Add patient IDs
    gene_expr_subset <- gene_expr_subset %>%
      left_join(sample_to_patient, by = "sample_id") %>%
      select(-sample_id)
    
    # Merge with clinical data
    merged_data <- gene_expr_subset %>%
      inner_join(clinical_data, by = "patient_id")
    
    cat("Merged data:", nrow(merged_data), "patients\n")
    
    if (nrow(merged_data) < 10) {
      cat("Too few samples for reliable multivariate analysis. Skipping.\n")
      next
    }
    
    # --- Standard Multivariate Cox ---
    
    # Create formula for multivariate Cox
    formula_str <- paste("Surv(overall_survival_time, overall_survival_status) ~", 
                        paste(significant_genes, collapse = " + "))
    formula <- as.formula(formula_str)
    
    # Fit multivariate Cox
    tryCatch({
      multi_cox <- coxph(formula, data = merged_data)
      
      # Summarize results
      multi_summary <- summary(multi_cox)
      
      # Extract coefficients
      multi_results <- data.frame(
        gene = rownames(multi_summary$coefficients),
        coef = multi_summary$coefficients[, "coef"],
        hazard_ratio = exp(multi_summary$coefficients[, "coef"]),
        ci_lower = exp(multi_summary$coefficients[, "coef"] - 
                     1.96 * multi_summary$coefficients[, "se(coef)"]),
        ci_upper = exp(multi_summary$coefficients[, "coef"] + 
                     1.96 * multi_summary$coefficients[, "se(coef)"]),
        p_value = multi_summary$coefficients[, "Pr(>|z|)"],
        significance = ifelse(multi_summary$coefficients[, "Pr(>|z|)"] < p_value_threshold, 
                            "Significant", "Not Significant")
      )
      
      # Adjust p-values
      multi_results$fdr <- p.adjust(multi_results$p_value, method = "BH")
      
      # Save results
      results_file <- file.path(results_dir, paste0("multivariate_cox_", cell_type, "_", comparison_name, ".csv"))
      write.csv(multi_results, results_file, row.names = FALSE)
      cat("Saved multivariate Cox results to", results_file, "\n")
      
      # Create forest plot for multivariate results
      forest_plot <- ggplot(multi_results, aes(x = reorder(gene, hazard_ratio), y = hazard_ratio)) +
        geom_point(aes(color = significance), size = 3) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = significance), width = 0.2) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "grey")) +
        coord_flip() +
        labs(
          title = paste("Multivariate Cox Regression -", comparison_name),
          x = "Gene",
          y = "Hazard Ratio (95% CI)",
          subtitle = paste("Cell type:", cell_type)
        ) +
        theme_bw()
      
      # Save plot
      forest_plot_file <- file.path(plot_dir, paste0("multivariate_forest_plot_", cell_type, "_", comparison_name, ".pdf"))
      ggsave(forest_plot_file, forest_plot, width = 10, height = 8)
      cat("Saved multivariate forest plot to", forest_plot_file, "\n")
      
    }, error = function(e) {
      cat("Error in multivariate Cox regression:", e$message, "\n")
    })
    
    # --- LASSO Cox Regression ---
    
    if (length(significant_genes) >= 3) {
      cat("Running LASSO Cox regression...\n")
      
      tryCatch({
        # Prepare data for glmnet
        x <- as.matrix(merged_data[, significant_genes])
        y <- Surv(merged_data$overall_survival_time, merged_data$overall_survival_status)
        
        # Fit LASSO model with cross-validation
        set.seed(123) # For reproducibility
        cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
        
        # Get optimal lambda
        optimal_lambda <- cv_fit$lambda.min
        cat("Optimal lambda:", optimal_lambda, "\n")
        
        # Fit model with optimal lambda
        lasso_fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = optimal_lambda)
        
        # Extract non-zero coefficients
        lasso_coefs <- coef(lasso_fit)
        nonzero_coefs <- lasso_coefs[which(lasso_coefs != 0), , drop = FALSE]
        
        if (length(nonzero_coefs) > 0) {
          # Convert to data frame
          lasso_results <- data.frame(
            gene = names(nonzero_coefs),
            coefficient = as.numeric(nonzero_coefs),
            hazard_ratio = exp(as.numeric(nonzero_coefs))
          )
          
          # Save results
          results_file <- file.path(results_dir, paste0("lasso_cox_", cell_type, "_", comparison_name, ".csv"))
          write.csv(lasso_results, results_file, row.names = FALSE)
          cat("Saved LASSO Cox results to", results_file, "\n")
          
          # Plot LASSO coefficients
          lasso_plot <- ggplot(lasso_results, aes(x = reorder(gene, coefficient), y = coefficient)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(
              title = paste("LASSO Cox Regression Coefficients -", comparison_name),
              x = "Gene",
              y = "Coefficient",
              subtitle = paste("Cell type:", cell_type, "| Lambda:", round(optimal_lambda, 4))
            ) +
            theme_bw()
          
          # Save plot
          lasso_plot_file <- file.path(plot_dir, paste0("lasso_coef_plot_", cell_type, "_", comparison_name, ".pdf"))
          ggsave(lasso_plot_file, lasso_plot, width = 10, height = 8)
          cat("Saved LASSO coefficient plot to", lasso_plot_file, "\n")
          
          # Create gene signature based on LASSO
          if (nrow(lasso_results) > 0) {
            cat("Creating gene signature based on LASSO coefficients...\n")
            
            # Calculate risk score for each patient
            risk_scores <- as.matrix(x[, lasso_results$gene, drop = FALSE]) %*% as.matrix(lasso_results$coefficient)
            
            # Add to data
            merged_data$risk_score <- as.numeric(risk_scores)
            
            # Dichotomize into high and low risk
            merged_data$risk_group <- ifelse(merged_data$risk_score > median(merged_data$risk_score), 
                                            "High Risk", "Low Risk")
            
            # Fit survival curve based on risk groups
            km_fit <- survfit(Surv(overall_survival_time, overall_survival_status) ~ risk_group, 
                             data = merged_data)
            
            # Plot K-M curve for risk groups
            risk_km_plot <- ggsurvplot(
              km_fit,
              data = merged_data,
              risk.table = TRUE,
              pval = TRUE,
              conf.int = TRUE,
              title = paste("Survival by LASSO Risk Score -", comparison_name),
              ggtheme = theme_bw(),
              risk.table.y.text.col = TRUE,
              risk.table.y.text = FALSE,
              palette = c("blue", "red")
            )
            
            # Save plot
            risk_km_file <- file.path(plot_dir, paste0("risk_score_km_", cell_type, "_", comparison_name, ".pdf"))
            pdf(risk_km_file, width = 10, height = 8)
            print(risk_km_plot)
            dev.off()
            cat("Saved risk score K-M plot to", risk_km_file, "\n")
            
            # Save risk scores
            risk_data <- merged_data %>%
              select(patient_id, risk_score, risk_group)
            
            risk_file <- file.path(results_dir, paste0("risk_scores_", cell_type, "_", comparison_name, ".csv"))
            write.csv(risk_data, risk_file, row.names = FALSE)
            cat("Saved risk scores to", risk_file, "\n")
          }
        } else {
          cat("No non-zero coefficients found in LASSO regression\n")
        }
        
      }, error = function(e) {
        cat("Error in LASSO Cox regression:", e$message, "\n")
      })
    } else {
      cat("Too few significant genes for LASSO regression. Skipping.\n")
    }
  }
}

# --- Proportional Hazards Assumption Testing ---

check_proportional_hazards <- function(cell_type, univariate_results) {
  cat("\n=== Checking Proportional Hazards Assumption for", cell_type, "===\n")
  
  if (is.null(univariate_results) || length(univariate_results) == 0) {
    cat("No univariate results available. Skipping proportional hazards check.\n")
    return(NULL)
  }
  
  # Process each comparison
  for (comparison_name in names(univariate_results)) {
    cat("Processing comparison:", comparison_name, "\n")
    
    # Get univariate results
    uni_results <- univariate_results[[comparison_name]]
    
    # Select top significant genes
    top_genes <- uni_results %>%
      filter(p_value < p_value_threshold) %>%
      head(top_n_genes) %>%
      pull(gene)
    
    if (length(top_genes) == 0) {
      cat("No significant genes for PH check in", comparison_name, "\n")
      next
    }
    
    cat("Testing proportional hazards assumption for", length(top_genes), "genes\n")
    
    # Subset TCGA expression data
    gene_expr_subset <- expression_data_final[top_genes, , drop = FALSE]
    gene_expr_subset <- t(gene_expr_subset)
    gene_expr_subset <- as.data.frame(gene_expr_subset) %>%
      rownames_to_column("sample_id")
    
    # Add patient IDs
    gene_expr_subset <- gene_expr_subset %>%
      left_join(sample_to_patient, by = "sample_id") %>%
      select(-sample_id)
    
    # Merge with clinical data
    merged_data <- gene_expr_subset %>%
      inner_join(clinical_data, by = "patient_id")
    
    # Initialize results
    ph_results <- data.frame()
    
    for (gene in top_genes) {
      tryCatch({
        # Create formula
        formula <- as.formula(paste("Surv(overall_survival_time, overall_survival_status) ~", gene))
        
        # Fit Cox model
        cox_model <- coxph(formula, data = merged_data)
        
        # Test proportional hazards
        ph_test <- cox.zph(cox_model)
        
        # Store results
        result <- data.frame(
          gene = gene,
          rho = ph_test$table[1, "rho"],
          chisq = ph_test$table[1, "chisq"],
          p_value = ph_test$table[1, "p"],
          proportional = ifelse(ph_test$table[1, "p"] > 0.05, "Yes", "No")
        )
        
        ph_results <- rbind(ph_results, result)
        
        # Plot Schoenfeld residuals for genes that violate the PH assumption
        if (ph_test$table[1, "p"] < 0.05) {
          cat("Creating Schoenfeld residual plot for gene", gene, "(violates PH assumption)\n")
          
          # Create plot
          residual_plot_file <- file.path(plot_dir, paste0("schoenfeld_residual_", cell_type, "_", 
                                                          comparison_name, "_", gene, ".pdf"))
          pdf(residual_plot_file, width = 8, height = 6)
          plot(ph_test)
          dev.off()
          cat("Saved Schoenfeld residual plot to", residual_plot_file, "\n")
        }
        
      }, error = function(e) {
        cat("Error testing PH assumption for gene", gene, ":", e$message, "\n")
      })
    }
    
    if (nrow(ph_results) > 0) {
      # Save results
      ph_file <- file.path(results_dir, paste0("proportional_hazards_", cell_type, "_", comparison_name, ".csv"))
      write.csv(ph_results, ph_file, row.names = FALSE)
      cat("Saved proportional hazards test results to", ph_file, "\n")
    }
  }
}

# --- Main Execution ---

# Run survival analysis for each cell type
for (cell_type in cell_types) {
  cat("\n\n=============================================================\n")
  cat("Starting survival analysis for", cell_type, "cell type\n")
  cat("=============================================================\n\n")
  
  # 1. Univariate Cox regression
  univariate_results <- run_univariate_cox(cell_type)
  
  if (!is.null(univariate_results) && length(univariate_results) > 0) {
    # 2. Multivariate Cox regression (incl. LASSO)
    run_multivariate_cox(cell_type, univariate_results)
    
    # 3. Check proportional hazards assumption
    check_proportional_hazards(cell_type, univariate_results)
  } else {
    cat("No univariate results available for", cell_type, ". Skipping further analysis.\n")
  }
}

cat("\nSurvival analysis with TCGA-OV data complete.\n")
cat("Results and plots saved in", output_dir, "\n") 