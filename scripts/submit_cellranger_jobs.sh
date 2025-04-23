#!/bin/bash

# Path to the metadata file
METADATA_FILE="/gpfs/commons/home/jameslee/HGSOC/metadata/GSE184880_metadata.csv"
OUTPUT_DIR="/gpfs/commons/home/jameslee/HGSOC/output/cellranger"
SCRIPT_DIR="/gpfs/commons/home/jameslee/HGSOC/scripts"
LOG_DIR="/gpfs/commons/home/jameslee/HGSOC/output/logs"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"
mkdir -p "$OUTPUT_DIR"

# Check if metadata file exists
if [ ! -f "$METADATA_FILE" ]; then
    echo "Error: Metadata file not found at $METADATA_FILE"
    exit 1
fi

# Function to submit a job for a specific run
submit_job() {
    local run_id=$1
    local sample_name=$2
    
    echo "Submitting job for run: $run_id (Sample: $sample_name)"
    
    # Create a job-specific script
    cat > "${SCRIPT_DIR}/job_${run_id}.sh" << EOL
#!/bin/bash
#SBATCH --job-name=cr_${run_id}
#SBATCH --output=${LOG_DIR}/${run_id}.out
#SBATCH --error=${LOG_DIR}/${run_id}.err
#SBATCH --time=48:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

# Path configurations
FASTQ_DIR="/gpfs/commons/home/jameslee/HGSOC/fastq/GSE184880"
OUTPUT_DIR="${OUTPUT_DIR}"

# Load necessary modules (adjust according to your system)
module load cellranger

# Create a temporary directory for this run's fastq files
temp_fastq_dir="\${OUTPUT_DIR}/tmp_${run_id}"
mkdir -p "\${temp_fastq_dir}"

# Create symlinks with CellRanger expected naming convention
ln -sf "\${FASTQ_DIR}/${run_id}_1.fastq" "\${temp_fastq_dir}/${sample_name}_S1_L001_R1_001.fastq"
ln -sf "\${FASTQ_DIR}/${run_id}_2.fastq" "\${temp_fastq_dir}/${sample_name}_S1_L001_R2_001.fastq"

# Run CellRanger count
cellranger count \\
    --id="${run_id}" \\
    --fastqs="\${temp_fastq_dir}" \\
    --sample="${sample_name}" \\
    --transcriptome=/path/to/reference/transcriptome \\
    --localcores=\$SLURM_CPUS_PER_TASK \\
    --localmem=90 \\
    --jobmode=local \\
    --disable-ui

# Move the output to the final directory
if [ -d "${run_id}/outs" ]; then
    mkdir -p "\${OUTPUT_DIR}/${run_id}"
    mv "${run_id}/outs" "\${OUTPUT_DIR}/${run_id}/"
    rm -rf "${run_id}"
fi

# Clean up temp directory
rm -rf "\${temp_fastq_dir}"

echo "Completed processing run: ${run_id}"
EOL

    # Make the job script executable
    chmod +x "${SCRIPT_DIR}/job_${run_id}.sh"
    
    # Submit the job
    sbatch "${SCRIPT_DIR}/job_${run_id}.sh"
}

# Process each run from the metadata file
# Skip the header line
tail -n +2 "$METADATA_FILE" | while IFS=',' read -r run_id sample_name rest_of_line; do
    submit_job "$run_id" "$sample_name"
done

echo "All jobs submitted." 