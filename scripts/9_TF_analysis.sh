#!/bin/bash
#SBATCH --job-name=TF_analysis
#SBATCH --output=TF_analysis_%j.out
#SBATCH --error=TF_analysis_%j.err
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load modules 
module load rstudio

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the integration script
Rscript 9_TF_analysis.R 

echo "TF job completed" 