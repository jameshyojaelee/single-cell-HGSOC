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

# Load required modules in correct order
module load GCCcore/12.3.0
module load OpenBLAS/0.3.23-GCC-12.3.0
module load GCC/12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0
module load FFTW/3.3.10-GCC-12.3.0
module load FlexiBLAS/3.3.1-GCC-12.3.0
module load ScaLAPACK/2.2.0-gompi-2023a-fb
module load foss/2023a
module load R/4.4.1
module load JAGS/4.3.2-foss-2023a

# Go to project directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the integration script
Rscript 11_infercnv_analysis.R 

echo "InferCNV job completed" 