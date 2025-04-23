#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --output=cellranger_%A_%a.out
#SBATCH --error=cellranger_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jameslee@nygenome.org

# Define GSE IDs
GSE_IDS=(
  "GSE193371"
  "GSE194105"
)

# Get current GSE ID
GSE_ID=${GSE_IDS[$SLURM_ARRAY_TASK_ID]}
echo "Processing $GSE_ID..."

# Set paths
CELLRANGER_PATH="/gpfs/commons/home/jameslee/cellranger-9.0.0"
FASTQ_DIR="/gpfs/commons/home/jameslee/HGSOC/fastq/${GSE_ID}"
OUTPUT_DIR="/gpfs/commons/home/jameslee/HGSOC/cellranger_output/${GSE_ID}"
REFERENCE="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2020-A" # Update this path to your reference

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Find all sample directories in the FASTQ directory
cd ${FASTQ_DIR}
SAMPLES=($(find . -name "*.fastq*" | grep -o "SRR[0-9]*" | sort | uniq))

# Check if samples were found
if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "No samples found in ${FASTQ_DIR}. Make sure to run fastq_download.sh first."
    exit 1
fi

echo "Found ${#SAMPLES[@]} samples to process."

# Function to run cellranger count for a sample
run_cellranger() {
    local sample=$1
    local output_dir="${OUTPUT_DIR}/${sample}"
    
    echo "Processing sample: ${sample}"
    echo "Fastq directory: ${FASTQ_DIR}"
    echo "Output directory: ${output_dir}"
    
    # Prepare fastq directory structure for cellranger
    # Create a directory structure cellranger expects
    TEMP_FASTQ_DIR="${OUTPUT_DIR}/fastq_temp/${sample}"
    mkdir -p ${TEMP_FASTQ_DIR}
    
    # Link fastq files to the expected format: sample_S1_L001_R1_001.fastq.gz
    ln -sf ${FASTQ_DIR}/${sample}_1.fastq.gz ${TEMP_FASTQ_DIR}/${sample}_S1_L001_R1_001.fastq.gz
    ln -sf ${FASTQ_DIR}/${sample}_2.fastq.gz ${TEMP_FASTQ_DIR}/${sample}_S1_L001_R2_001.fastq.gz
    
    # Run cellranger
    ${CELLRANGER_PATH}/cellranger count \
        --id=${sample} \
        --transcriptome=${REFERENCE} \
        --fastqs=${TEMP_FASTQ_DIR} \
        --sample=${sample} \
        --localcores=16 \
        --localmem=120         
    # Move results to final output directory
    mkdir -p ${output_dir}
    mv ${sample} ${output_dir}/
    
    echo "Completed processing for sample: ${sample}"
}

# Process each sample
for sample in "${SAMPLES[@]}"; do
    run_cellranger ${sample}
done

echo "All samples for ${GSE_ID} have been processed." 