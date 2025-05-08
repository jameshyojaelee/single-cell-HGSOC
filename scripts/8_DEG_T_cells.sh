#!/bin/bash
#SBATCH --job-name=DEG_T_cells
#SBATCH --output=DEG_T_cells_%j.out
#SBATCH --error=DEG_T_cells_%j.err
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=210G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load modules 
module load rstudio

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the T cells DEG script
Rscript 8_DEG_T_cells.R

echo "DEG_T_cells job completed"