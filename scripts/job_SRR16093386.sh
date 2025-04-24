#!/bin/bash
#SBATCH --job-name=cr_SRR16093386
#SBATCH --output=/gpfs/commons/home/jameslee/HGSOC/output/logs/SRR16093386.out
#SBATCH --error=/gpfs/commons/home/jameslee/HGSOC/output/logs/SRR16093386.err
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --partition=cpu
#SBATCH --mail-type=END
#SBATCH --mail-user=jameslee@nygenome.org

# Path configurations
FASTQ_DIR="/gpfs/commons/home/jameslee/HGSOC/fastq/GSE184880"
OUTPUT_DIR="/gpfs/commons/home/jameslee/HGSOC/output/cellranger"

# Create a temporary directory for this run's fastq files
temp_fastq_dir="${OUTPUT_DIR}/tmp_SRR16093386"
mkdir -p "${temp_fastq_dir}"

# Create symlinks with CellRanger expected naming convention
ln -sf "${FASTQ_DIR}/SRR16093386_1.fastq" "${temp_fastq_dir}/GSM5599231_S1_L001_R1_001.fastq"
ln -sf "${FASTQ_DIR}/SRR16093386_2.fastq" "${temp_fastq_dir}/GSM5599231_S1_L001_R2_001.fastq"

# Run CellRanger count
/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger count \
    --id="SRR16093386" \
    --fastqs="${temp_fastq_dir}" \
    --sample="GSM5599231" \
    --transcriptome=/gpfs/commons/home/acorman/sanjanalab/ASD_inducible_screen/Pilot_QCMiseq/GRCh38/refdata-gex-GRCh38-2024-A \
    --create-bam=true \
    --localcores=$SLURM_CPUS_PER_TASK \
    --localmem=60 \
    --jobmode=local \
    --disable-ui

# Move the output to the final directory
if [ -d "SRR16093386/outs" ]; then
    mkdir -p "${OUTPUT_DIR}/SRR16093386"
    mv "SRR16093386/outs" "${OUTPUT_DIR}/SRR16093386/"
    rm -rf "SRR16093386"
fi

# Clean up temp directory
rm -rf "${temp_fastq_dir}"

echo "Completed processing run: SRR16093386"
