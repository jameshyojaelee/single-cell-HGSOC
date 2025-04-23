#!/bin/bash
#SBATCH --job-name=fastq_download2
#SBATCH --output=fastq_download2_%A_%a.out
#SBATCH --error=fastq_download2_%A_%a.err
#SBATCH --time=90:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Load required environment and modules
source /gpfs/commons/home/jameslee/miniconda3/etc/profile.d/conda.sh
conda activate RNA
module load SRA-Toolkit/3.2.0-gompi-2023b

# Path to the Python script
SCRIPT="/gpfs/commons/home/jameslee/HGSOC/fastq_download2.py"

# Make sure script is executable
chmod +x ${SCRIPT}

# Array tasks:
# 0: GSE184880 (full download)
# 1-5: GSE180661 (chunked download, 5 chunks)

# Define number of chunks for GSE180661
NUM_CHUNKS=5

# Print job info
echo "Running job ${SLURM_ARRAY_TASK_ID}"

if [ ${SLURM_ARRAY_TASK_ID} -eq 0 ]; then
    # Process GSE184880 (no chunking)
    echo "Processing GSE184880 (full dataset)"
    ${SCRIPT} GSE184880
else
    # Process GSE180661 in chunks
    # Calculate chunk index (0-based)
    CHUNK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    echo "Processing GSE180661 (chunk ${CHUNK_INDEX} of ${NUM_CHUNKS})"
    ${SCRIPT} GSE180661 --chunk ${CHUNK_INDEX} --total-chunks ${NUM_CHUNKS}
fi 