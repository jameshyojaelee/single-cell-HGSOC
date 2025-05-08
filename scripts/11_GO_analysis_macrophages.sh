#!/bin/bash
#SBATCH --job-name=GO_macrophages
#SBATCH --output=GO_macrophages_%j.out
#SBATCH --error=GO_macrophages_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load modules 
module load rstudio

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the macrophages GO analysis script
Rscript 11_GO_analysis_macrophages.R

echo "GO_macrophages job completed" 