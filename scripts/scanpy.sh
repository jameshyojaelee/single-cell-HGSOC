#!/bin/bash
#SBATCH --job-name=scanpy_process
#SBATCH --output=scanpy_process_%j.out
#SBATCH --error=scanpy_process_%j.err
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=196G
#SBATCH --partition=cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=jameslee@nygenome.org

# Load required modules
module purge
source $(conda info --base)/etc/profile.d/conda.sh
conda activate scanpy

# Set working directory
cd /gpfs/commons/home/jameslee/HGSOC/scripts

# Run the Python script
python scanpy_integration.py

echo "Scanpy processing job completed at $(date)"