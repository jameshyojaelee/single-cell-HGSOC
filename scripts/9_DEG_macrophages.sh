#!/bin/bash
#SBATCH --job-name=DEG_macrophages
#SBATCH --output=DEG_macrophages_%j.out
#SBATCH --error=DEG_macrophages_%j.err
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

# Run the macrophages DEG script
Rscript 9_DEG_macrophages.R

echo "DEG_macrophages job completed" 