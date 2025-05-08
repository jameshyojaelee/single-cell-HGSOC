#!/bin/bash
#SBATCH --job-name=TCGA_Survival
#SBATCH --output=TCGA_Survival_%j.out
#SBATCH --error=TCGA_Survival_%j.err
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load modules 
module load rstudio

# Set error handling to capture and display errors more clearly
set -e

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Check and install required packages if needed
Rscript -e '
pkg_check <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    if (pkg %in% c("TCGAbiolinks", "SummarizedExperiment", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "fgsea", "ComplexHeatmap")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg, repos="http://cran.us.r-project.org")
    }
  }
}

# List of required packages
pkgs <- c("survival", "survminer", "dplyr", "ggplot2", "readr", "tidyr", "TCGAbiolinks", "glmnet", 
          "stringr", "ComplexHeatmap", "RColorBrewer", "SummarizedExperiment", "fgsea")

# Check and install if needed
for (pkg in pkgs) pkg_check(pkg)
'

# Run the TCGA survival analysis script
echo "Starting TCGA survival analysis at $(date)"
Rscript 12_Survival_Analysis_TCGA.R

if [ $? -eq 0 ]; then
    echo "TCGA survival analysis job completed successfully at $(date)"
else
    echo "ERROR: TCGA survival analysis job failed at $(date)"
    echo "Check the error log for details"
    exit 1
fi