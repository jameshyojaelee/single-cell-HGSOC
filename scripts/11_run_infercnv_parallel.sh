#!/bin/bash
#SBATCH --job-name=InferCNV_Array
#SBATCH --output=infercnv_array_%A_%a.out  # %A = job ID, %a = array task ID
#SBATCH --error=infercnv_array_%A_%a.err
#SBATCH --time=24:00:00  # Reduced time estimate per task
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # Matches R script's num_threads
#SBATCH --mem=128G      
#SBATCH --array=0-2 
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS
#SBATCH --mail-user=jameslee@nygenome.org

# =====================================================
# Run InferCNV using a SLURM Job Array
# =====================================================

set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.

# Define the path to the R script runner
R_SCRIPT="HGSOC/scripts/11_infercnv_stage_runner.R"

# Define the tumor stages corresponding to array indices
# Index 0 -> IC2, Index 1 -> IIB, etc.
STAGES=("IIB" "IIIB" "IIIC")

# Get the stage for the current array task
if [ -z "${SLURM_ARRAY_TASK_ID+x}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set. This script must be run with sbatch." >&2
    exit 1
fi
TASK_ID=${SLURM_ARRAY_TASK_ID}
CURRENT_STAGE=${STAGES[$TASK_ID]}

# Ensure the R script exists
if [ ! -f "${R_SCRIPT}" ]; then
    echo "Error: R script not found at ${R_SCRIPT}" >&2
    exit 1
fi

echo "Starting InferCNV task ${TASK_ID} for stage: ${CURRENT_STAGE}"

# Create necessary log directory for SLURM output files
mkdir -p HGSOC/logs/slurm

# Execute the R script for the assigned stage
# R script output (cat messages) will go to the SLURM .out file
# R script errors/warnings will go to the SLURM .err file
Rscript "${R_SCRIPT}" "${CURRENT_STAGE}"

exit_code=$?

echo "Finished InferCNV task ${TASK_ID} for stage: ${CURRENT_STAGE} with exit code ${exit_code}"

exit ${exit_code} 