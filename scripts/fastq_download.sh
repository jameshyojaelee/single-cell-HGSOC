#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=fastq_download_%A_%a.out
#SBATCH --error=fastq_download_%A_%a.err
#SBATCH --time=90:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load required environment and modules
source /gpfs/commons/home/jameslee/miniconda3/etc/profile.d/conda.sh
conda activate RNA
module load SRA-Toolkit/3.2.0-gompi-2023b

# Define GSE IDs for each array task
GSE_IDS=(
  "GSE193371"
  "GSE194105"
)

# Path to the Python script
SCRIPT="/gpfs/commons/home/jameslee/HGSOC/scripts/fastq_download.py"

# Make sure script is executable
chmod +x ${SCRIPT}

# Print job info
echo "Running job ${SLURM_ARRAY_TASK_ID} for ${GSE_IDS[$SLURM_ARRAY_TASK_ID]}"

# Run the Python script with the appropriate GSE ID
${SCRIPT} ${GSE_IDS[$SLURM_ARRAY_TASK_ID]}