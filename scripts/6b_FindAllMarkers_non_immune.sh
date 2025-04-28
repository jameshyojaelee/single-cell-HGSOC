#!/bin/bash
#SBATCH --job-name=FindAllMarkers_non_immune
#SBATCH --output=FindAllMarkers_non_immune_%j.out
#SBATCH --error=FindAllMarkers_non_immune_%j.err
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
Rscript 6b_FindAllMarkers_non_immune.R

echo "FindAllMarkers_non_immune job completed" 