#!/bin/bash
#SBATCH --job-name=seurat_full
#SBATCH --output=seurat_full_%j.out
#SBATCH --error=seurat_full_%j.err
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load required modules
module purge
module load rstudio

# Set working directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the combined R script
Rscript 2_process_and_integrate.R

echo "Combined pipeline job completed at $(date)" 