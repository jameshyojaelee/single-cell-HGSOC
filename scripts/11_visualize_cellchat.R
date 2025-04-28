# ======================================================
# Visualize CellChat Results from Saved Objects
# ======================================================

# Load required libraries
library(CellChat, lib="/gpfs/commons/home/jameslee/R/x86_64-pc-linux-gnu-library/4.4/")
library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)
library(stringr)

# --- Configuration ---

# Set working directory and paths
base_dir <- "/gpfs/commons/home/jameslee/HGSOC"
cellchat_output_dir <- file.path(base_dir, "output/cellchat")

# Check if base output directory exists
if (!dir.exists(cellchat_output_dir)) {
  stop(paste("CellChat output directory not found:", cellchat_output_dir, 
             "\nPlease run script 10_cellchat_analysis.R first."))
}

# --- Visualize Single Stage Results ---

cat("Processing individual stage CellChat results...\n")

# Find the stage-specific RDS files (adjust pattern if needed)
rds_files <- list.files(path = cellchat_output_dir, 
                        pattern = "cellchat_grouped_.*\\.rds", 
                        full.names = TRUE,
                        recursive = TRUE) # Search in subdirectories too

# Filter out the merged file if it's caught by the pattern
rds_files <- rds_files[!grepl("comparison_grouped", rds_files)]

# Initialize a list to store loaded CellChat objects for later combined plotting
object.list <- list()

if (length(rds_files) == 0) {
  cat("No stage-specific CellChat RDS files found matching the pattern.\n")
} else {
  cat(paste("Found", length(rds_files), "stage-specific RDS files.\n"))
  
  for (rds_path in rds_files) {
    stage_name <- str_extract(basename(rds_path), "(?<=cellchat_grouped_)(.*)(?=\\.rds)")
    stage_output_dir <- file.path(cellchat_output_dir, stage_name)
    
    cat(paste("\nProcessing stage:", stage_name, "from", rds_path, "\n"))
    
    if (!dir.exists(stage_output_dir)) {
       cat(paste("  Warning: Stage output directory not found, creating:", stage_output_dir, "\n"))
       dir.create(stage_output_dir, recursive = TRUE)
    }
    
    tryCatch({
      # Load the CellChat object
      cat("  Loading CellChat object...")
      cellchat <- readRDS(rds_path)
      cat(" Done.\n")
      
      # Ensure necessary slots exist (simple check)
      if (is.null(cellchat@net) || is.null(cellchat@net$weight)) {
          stop("Required slots (@net or @net$weight) not found in the loaded object.")
      }
      if (is.null(cellchat@idents)) {
          stop("Required slot @idents not found.")
      }
      
      # Calculate groupSize if not stored or needs recalculation (based on idents)
      # It should ideally be correct from the analysis script, but recalculate for safety
      groupSize <- as.numeric(table(cellchat@idents))
      
      # Add the loaded object to the list for later use
      object.list[[stage_name]] <- cellchat
      
      # Define plot file path
      plot_file_path <- file.path(stage_output_dir, paste0("plots_grouped_", stage_name, ".pdf"))
      pdf(plot_file_path, width = 10, height = 8)
      cat(paste("  Generating plots and saving to:", plot_file_path, "\n"))
      
      # Interaction weights/strength heatmap
      tryCatch({
          print(CellChat::netVisual_heatmap(cellchat, measure = "weight", 
                                            title = paste("Interaction Weights (Grouped) -", stage_name)))
      }, error=function(e){cat(paste("    Error generating heatmap plot:", e$message, "\n"))})
      
      # Circle plot
      tryCatch({
          print(netVisual_circle(cellchat@net$weight, 
                                 vertex.weight = groupSize, weight.scale = T, 
                                 label.edge= F, title.name = paste("Interaction Weights/Number (Grouped) -", stage_name)))
      }, error=function(e){cat(paste("    Error generating circle plot:", e$message, "\n"))})
      
      # Add other desired plots here (e.g., pathway analysis, L-R pairs)
      # Example: Visualize top pathways
      # tryCatch({
      #    cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional") # Example step if needed
      #    cellchat <- netEmbedding(cellchat, type = "functional")
      #    cellchat <- netClustering(cellchat, type = "functional")
      #    print(plotGeneExpression(cellchat, features = c("TGFB1")))
      # }, error=function(e){cat(paste("    Error generating extra plots:", e$message, "\n"))})
      
      dev.off()
      cat(paste("  Finished plots for stage", stage_name, "\n"))
      
    }, error = function(e) {
      cat(paste("  Error processing file", rds_path, ":", e$message, "\n"))
    })
  }
}

# --- Generate Combined Signaling Role Scatter Plots ---

cat("\nGenerating combined signaling role scatter plots across stages...\n")

if (length(object.list) < 2) { # Need at least two objects to compare roles meaningfully, adjust if needed
  cat("Skipping signaling role scatter plots: Need at least 2 successfully processed stages.\n")
} else {
  tryCatch({
    # Calculate link numbers and min/max weights across all datasets
    num.link <- sapply(object.list, function(x) {
      if (!is.null(x@net$count)) {
        rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
      } else {
        numeric(0) # Return empty if counts are missing
      }
    })
    # Ensure num.link is a vector or matrix, not a list of empty numerics if errors occurred
    num.link_vector <- unlist(num.link)
    
    if (length(num.link_vector) > 0) {
      weight.MinMax <- c(min(num.link_vector), max(num.link_vector))
      cat(paste("  Calculated MinMax weight for scatter plots:", paste(weight.MinMax, collapse=" - "), "\n"))
      
      gg <- list()
      for (i in 1:length(object.list)) {
        stage_title <- names(object.list)[i]
        cat(paste("  Generating scatter plot for:", stage_title, "\n"))
        tryCatch({
          # Compute centrality scores before plotting
          cat("    Computing centrality...")
          current_object <- object.list[[i]]
          current_object <- netAnalysis_computeCentrality(current_object, slot.name = "netP") 
          cat(" Done.\n")
          # Update the object in the list with the computed centrality
          object.list[[i]] <- current_object 
          
          # Now generate the plot using the updated object
          gg[[i]] <- netAnalysis_signalingRole_scatter(current_object, 
                                                       title = stage_title, 
                                                       weight.MinMax = weight.MinMax)
        }, error = function(e) {
          cat(paste("    Error generating scatter plot for", stage_title, ":", e$message, "\n"))
          gg[[i]] <- NULL # Add placeholder if error occurs
        })
      }
      
      # Filter out NULL plots from gg list if errors occurred
      gg <- gg[!sapply(gg, is.null)]
      
      if (length(gg) > 0) {
        # Save the combined plot
        scatter_plot_path <- file.path(cellchat_output_dir, "plot_signaling_role_scatter_all_stages.pdf")
        pdf(scatter_plot_path, width = max(10, 4 * length(gg)), height = 10)
        cat(paste("  Saving combined scatter plots to:", scatter_plot_path, "\n"))
        print(patchwork::wrap_plots(plots = gg))
        dev.off()
      } else {
         cat("  Skipping scatter plot PDF: No scatter plots were successfully generated.\n")
      }
    } else {
      cat("  Skipping scatter plots: Could not calculate link numbers (perhaps @net$count missing?).\n")
    }
    
  }, error = function(e) {
    cat(paste("  Error during combined signaling role scatter plot generation:", e$message, "\n"))
  })
}

# --- Visualize Comparative Analysis Results ---

cat("\nProcessing comparative analysis results...\n")
comparison_dir <- file.path(cellchat_output_dir, "comparison_grouped")
merged_rds_path <- file.path(comparison_dir, "cellchat_merged_grouped.rds")

# Initialize merged object variable
cellchat_merged <- NULL

# Check if merged file exists, if not, create it
if (!file.exists(merged_rds_path)) {
  cat(paste("Merged CellChat object not found at:", merged_rds_path, "\n"))
  
  # Check if we have loaded individual objects in object.list
  if (length(object.list) >= 2) {
    cat("Attempting to create merged object from loaded individual stages...\n")
    tryCatch({
      if (!dir.exists(comparison_dir)) {
        dir.create(comparison_dir, recursive = TRUE)
      }
      cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
      cat("  Saving newly created merged object...")
      saveRDS(cellchat_merged, file = merged_rds_path)
      cat(" Done.\n")
    }, error = function(e) {
      cat(paste("  Error creating or saving merged object:", e$message, "\n"))
      cellchat_merged <- NULL # Ensure it remains NULL if creation failed
    })
  } else {
    cat("Skipping merge: Fewer than 2 individual stage objects were successfully loaded earlier.\n")
  }
  
} else {
  cat(paste("Found existing merged CellChat object:", merged_rds_path, "\n"))
  # Try loading the existing file
  tryCatch({
    cat("  Loading existing merged CellChat object...")
    cellchat_merged <- readRDS(merged_rds_path)
    cat(" Done.\n")
  }, error = function(e) {
     cat(paste("  Error loading existing merged file", merged_rds_path, ":", e$message, "\n"))
     cellchat_merged <- NULL
  })
}

# Proceed only if cellchat_merged is not NULL (i.e., loaded or created successfully)
if (!is.null(cellchat_merged)) {
  # Check if there are multiple datasets to compare
  if (length(levels(cellchat_merged@meta$datasets)) < 2) {
    cat("  Skipping comparison plots: Merged object does not contain multiple datasets/stages to compare.\n")
  } else {
    # Define comparison plot file path
    comparison_plot_path <- file.path(comparison_dir, "plots_comparison_grouped.pdf")
    pdf(comparison_plot_path, width = 12, height = 10)
    cat(paste("  Generating comparison plots and saving to:", comparison_plot_path, "\n"))
    
    # Compare total interactions
    tryCatch({
        group_names <- levels(cellchat_merged@meta$datasets)
        print(compareInteractions(cellchat_merged, show.legend = F, group = group_names,
                                  title = "Comparison: Total Interactions (Grouped)"))
    }, error=function(e){cat(paste("    Error generating compareInteractions plot:", e$message, "\n"))})
    
    # Compare interactions heatmap
    tryCatch({
         print(CellChat::netVisual_heatmap(cellchat_merged, measure = "weight", 
                                           title = "Comparison: Interaction Weights (Grouped)"))
    }, error=function(e){cat(paste("    Error generating comparison heatmap:", e$message, "\n"))})
    
    # Differential interaction plot (comparing first two groups)
    if(length(group_names) >= 2) {
        comparison_groups <- c(1, 2) # Define which groups to compare (e.g., first vs second)
        cat(paste("  Generating differential interaction plot comparing groups:", 
                  group_names[comparison_groups[1]], "vs", group_names[comparison_groups[2]], "\n"))
        tryCatch({
            # Remove slot.name argument
            print(netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", 
                                          comparison = comparison_groups))
        }, error=function(e){cat(paste("    Error generating differential interaction plot:", e$message, "\n"))})
    } else {
         cat("  Skipping differential interaction plot: Need at least two groups to compare.\n")
    }
    
    # --- Differential Pathway Analysis --- 
    cat("  Generating differential pathway analysis plots (rankNet)...")
    tryCatch({
        # Weight-based ranking
        print(rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE, 
                      title = "Differential Pathway Weight Ranking", measure = "weight"))
        # Count-based ranking (optional)
        # print(rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE, 
        #              title = "Differential Pathway Count Ranking", measure = "count"))
    }, error=function(e){cat(paste("    Error generating rankNet plots:", e$message, "\n"))})
    
    # Add other comparison plots as desired (e.g., differential network analysis)
    # Refer to CellChat comparison vignettes for examples
    
    dev.off()
    cat(paste("  Finished comparison plots.\n"))
  }
}

cat("\nCellChat visualization script complete.\n") 