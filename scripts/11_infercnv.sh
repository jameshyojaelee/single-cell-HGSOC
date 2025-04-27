#!/bin/bash
#SBATCH --job-name=infercnv_analysis
#SBATCH --output=infercnv_analysis_%j.out
#SBATCH --error=infercnv_analysis_%j.err
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Purge all modules first
module purge

# Load modules in correct order
module load rstudio
module load JAGS/4.3.2-foss-2023a

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the integration script
Rscript 11_infercnv_analysis.R 

echo "InferCNV job completed" 