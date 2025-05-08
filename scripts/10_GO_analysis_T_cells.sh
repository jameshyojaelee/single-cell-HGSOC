#!/bin/bash
#SBATCH --job-name=GO_T_cells
#SBATCH --output=GO_T_cells_%j.out
#SBATCH --error=GO_T_cells_%j.err
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

# Run the T cells GO analysis script
Rscript 10_GO_analysis_T_cells.R

echo "GO_T_cells job completed" 