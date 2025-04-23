#!/gpfs/commons/home/jameslee/miniconda3/envs/RNA/bin/python
import sys
import subprocess
import os
import argparse

def download_fastq(srr_id, output_dir):
    print(f"Downloading fastq for {srr_id}...")
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"cd {output_dir} && fastq-dump --split-files {srr_id}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Successfully downloaded {srr_id}")
    except Exception as e:
        print(f"Error downloading {srr_id}: {e}")

def query_and_download(gse_id):
    try:
        from pysradb import SRAweb
    except ImportError:
        sys.exit("pysradb not found. Install via: pip install pysradb")

    output_dir = f"/gpfs/commons/home/jameslee/HGSOC/fastq/{gse_id}"
    
    db = SRAweb()
    print(f"Querying SRA for {gse_id}...")
    srp_df = db.gse_to_srp(gse_id)
    if srp_df.empty:
        print(f"No study accessions found for {gse_id}")
    else:
        for srp in srp_df['study_accession']:
            metadata = db.sra_metadata(srp, detailed=False)
            if metadata.empty:
                print(f"No run accessions found for study {srp}")
            else:
                for srr in metadata['run_accession']:
                    download_fastq(srr, output_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download FASTQ files from SRA for a GSE ID')
    parser.add_argument('gse_id', help='GSE ID to download')
    
    args = parser.parse_args()
    query_and_download(args.gse_id) 