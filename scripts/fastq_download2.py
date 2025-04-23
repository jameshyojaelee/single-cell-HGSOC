#!/gpfs/commons/home/jameslee/miniconda3/envs/RNA/bin/python
import sys
import subprocess
import os
import argparse
import math

def download_fastq(srr_id, output_dir):
    print(f"Downloading fastq for {srr_id}...")
    os.makedirs(output_dir, exist_ok=True)
    cmd = f"cd {output_dir} && fastq-dump --split-files {srr_id}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Successfully downloaded {srr_id}")
    except Exception as e:
        print(f"Error downloading {srr_id}: {e}")

def query_and_download(gse_id, chunk_index=None, total_chunks=None):
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
        return
    
    # Collect all SRR accessions first
    all_srr_ids = []
    for srp in srp_df['study_accession']:
        metadata = db.sra_metadata(srp, detailed=False)
        if metadata.empty:
            print(f"No run accessions found for study {srp}")
        else:
            all_srr_ids.extend(metadata['run_accession'].tolist())
    
    # Save all SRR IDs to a file for reference
    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/all_srr_ids.txt", "w") as f:
        for srr in all_srr_ids:
            f.write(f"{srr}\n")
    
    print(f"Found {len(all_srr_ids)} SRR accessions for {gse_id}")
    
    # If chunking is specified, only process the assigned chunk
    if chunk_index is not None and total_chunks is not None:
        chunk_size = math.ceil(len(all_srr_ids) / total_chunks)
        start_idx = chunk_index * chunk_size
        end_idx = min(start_idx + chunk_size, len(all_srr_ids))
        
        chunk_srr_ids = all_srr_ids[start_idx:end_idx]
        print(f"Processing chunk {chunk_index+1} of {total_chunks}: {len(chunk_srr_ids)} SRR accessions ({start_idx}-{end_idx-1})")
        
        # Download only the SRRs in this chunk
        for srr in chunk_srr_ids:
            download_fastq(srr, output_dir)
    else:
        # Download all SRRs
        for srr in all_srr_ids:
            download_fastq(srr, output_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download FASTQ files from SRA for a GSE ID')
    parser.add_argument('gse_id', choices=['GSE184880', 'GSE180661'], 
                       help='GSE ID to download (GSE184880 or GSE180661)')
    parser.add_argument('--chunk', type=int, help='Chunk index (0-based) to process')
    parser.add_argument('--total-chunks', type=int, help='Total number of chunks')
    
    args = parser.parse_args()
    
    # Both or none of chunk arguments must be provided
    if (args.chunk is not None) != (args.total_chunks is not None):
        parser.error("--chunk and --total-chunks must be provided together")
    
    # Validate chunk index if provided
    if args.chunk is not None and args.chunk < 0:
        parser.error("Chunk index must be non-negative")
    
    if args.total_chunks is not None and args.total_chunks <= 0:
        parser.error("Total chunks must be positive")
    
    query_and_download(args.gse_id, args.chunk, args.total_chunks) 