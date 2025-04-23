# CellRanger Count Pipeline for GSE184880

This directory contains scripts to run CellRanger count on the GSE184880 scRNA-seq dataset on an HPC cluster.

## Files

- `submit_cellranger_jobs.sh`: Script to submit SLURM jobs for parallel processing on HPC
- `job_*.sh`: Individual job scripts generated for each run (created when `submit_cellranger_jobs.sh` is executed)

## Prerequisites

1. CellRanger must be available via the module system on the HPC
2. Reference transcriptome must be available
3. Metadata file (`GSE184880_metadata.csv`) must be in the correct location
4. FASTQ files must be in the expected location (`fastq/GSE184880/`)

## Usage

To submit jobs to the SLURM scheduler:

```bash
bash submit_cellranger_jobs.sh
```

This will create and submit individual job scripts for each run, allowing parallel processing.

## Configuration

Before running the scripts, you should modify the following parameters:

1. Update the path to the reference transcriptome:
   ```bash
   --transcriptome=/path/to/reference/transcriptome
   ```

2. Adjust the resource requirements if needed:
   ```bash
   #SBATCH --time=48:00:00
   #SBATCH --mem=96G
   #SBATCH --cpus-per-task=16
   ```

3. Make sure the module name for CellRanger is correct for your system:
   ```bash
   module load cellranger
   ```

## Output

The output will be stored in the `output/cellranger` directory, with each run having its own subdirectory containing the CellRanger output files.

## Monitoring

The log files for each job will be stored in the `output/logs` directory. You can check the status of your jobs using:

```bash
squeue -u $USER
```

And check the output logs using:

```bash
tail -f output/logs/SRR*.out
``` 