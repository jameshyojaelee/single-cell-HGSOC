#!/bin/bash
#SBATCH --job-name=seurat_integrate
#SBATCH --output=seurat_integrate_%j.out
#SBATCH --error=seurat_integrate_%j.err
#SBATCH --time=60:00:00
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
Rscript 3_integrate_samples.R

echo "Integration job completed at $(date)"
