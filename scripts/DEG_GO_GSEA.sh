#!/bin/bash
#SBATCH --job-name=DEG_GO_GSEA
#SBATCH --output=DEG_GO_GSEA_%j.out
#SBATCH --error=DEG_GO_GSEA_%j.err
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
Rscript 6_DEG_non_immune.R
Rscript 7_GO_analysis.R
Rscript 8_GSEA_analysis.R

echo "DEG, GO, and GSEA job completed" 