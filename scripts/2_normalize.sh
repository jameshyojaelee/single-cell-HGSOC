#!/bin/bash
#SBATCH --job-name=seurat_norm
#SBATCH --output=seurat_norm_%j.out
#SBATCH --error=seurat_norm_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=196G
#SBATCH --partition=cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=jameslee@nygenome.org

# Load required modules
module purge
module load rstudio

# Set working directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the R script
Rscript 2_normalize.R

echo "Job completed at $(date)"