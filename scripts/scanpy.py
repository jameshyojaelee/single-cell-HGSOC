#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import defaultdict

# Set random seed for reproducibility
np.random.seed(42)

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, figsize=(8, 8))

# Define paths
BASE_DIR = "/gpfs/commons/home/jameslee/HGSOC"
CELLRANGER_DIR = os.path.join(BASE_DIR, "output/cellranger")
METADATA_PATH = os.path.join(BASE_DIR, "metadata/GSE184880_metadata.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output/scanpy")

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load metadata
metadata = pd.read_csv(METADATA_PATH)
print(f"Loaded metadata with {len(metadata)} entries")

# Group SRR runs by GSM sample name
sample_groups = defaultdict(list)
for _, row in metadata.iterrows():
    sample_groups[row['Sample Name']].append(row['Run'])

print(f"Found {len(sample_groups)} unique samples")

# Process each sample
all_samples = []
sample_adatas = {}

for sample_name, run_ids in sample_groups.items():
    print(f"\nProcessing {sample_name} with {len(run_ids)} sequencing runs")
    
    # List to store run-level anndata objects
    run_adatas = []
    
    # Process each run
    for run_id in run_ids:
        h5_path = os.path.join(CELLRANGER_DIR, run_id, "outs", "filtered_feature_bc_matrix.h5")
        
        if not os.path.exists(h5_path):
            print(f"Warning: CellRanger output for {run_id} not found at {h5_path}")
            continue
        
        # Load data
        print(f"Loading {run_id}...")
        adata = sc.read_10x_h5(h5_path)
        
        # Add metadata
        run_meta = metadata[metadata['Run'] == run_id].iloc[0]
        adata.obs['run_id'] = run_id
        adata.obs['sample_name'] = sample_name
        adata.obs['age'] = run_meta['age']
        adata.obs['pathology'] = run_meta['pathology'] 
        adata.obs['tissue_type'] = run_meta['tissue_type']
        adata.obs['tumor_stage'] = run_meta['tumor_stage']
        
        # Add unique barcode suffix to avoid conflicts when merging
        adata.obs_names = [f"{bc}_{run_id}" for bc in adata.obs_names]
        
        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        
        # Annotate mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        
        run_adatas.append(adata)
    
    if not run_adatas:
        print(f"No data found for {sample_name}, skipping")
        continue
    
    # Merge all runs for this sample if there are multiple runs
    if len(run_adatas) > 1:
        adata_combined = run_adatas[0].concatenate(run_adatas[1:], batch_key='sequencing_run')
        print(f"Merged {len(run_adatas)} runs for {sample_name}, with {adata_combined.n_obs} cells")
    else:
        adata_combined = run_adatas[0]
        adata_combined.obs['sequencing_run'] = '0'
        print(f"Single run for {sample_name} with {adata_combined.n_obs} cells")
    
    # QC filtering
    sc.pl.violin(adata_combined, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                jitter=0.4, multi_panel=True, save=f"_{sample_name}_qc_metrics.pdf")
    
    # Filter cells based on QC metrics
    adata_filtered = adata_combined[adata_combined.obs.n_genes_by_counts > 200, :]
    adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 40, :]
    
    print(f"After filtering: {adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes")
    
    # Normalize data
    sc.pp.normalize_total(adata_filtered, target_sum=1e4)
    sc.pp.log1p(adata_filtered)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Plot variable genes
    sc.pl.highly_variable_genes(adata_filtered, save=f"_{sample_name}_hvg.pdf")
    
    # Scale data
    sc.pp.scale(adata_filtered, max_value=10)
    
    # Run PCA
    sc.tl.pca(adata_filtered, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata_filtered, n_pcs=50, save=f"_{sample_name}_pca_variance.pdf")
    
    # Compute neighborhood graph and UMAP
    sc.pp.neighbors(adata_filtered, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata_filtered)
    
    # Clustering
    sc.tl.leiden(adata_filtered, resolution=0.5)
    
    # Plot UMAP
    sc.pl.umap(adata_filtered, color=['leiden', 'sequencing_run', 'tissue_type', 'tumor_stage'], 
              wspace=0.4, save=f"_{sample_name}_umap.pdf")
    
    # Known cell type markers
    markers = {
        'Epithelial': ['EPCAM', 'KRT8', 'KRT18'],
        'Immune': ['PTPRC', 'CD3E', 'CD4', 'CD8A', 'CD14', 'CD68'],
        'Fibroblast': ['COL1A1', 'DCN', 'LUM']
    }
    
    # Plot markers if they exist in the dataset
    for cell_type, marker_genes in markers.items():
        existing_markers = [m for m in marker_genes if m in adata_filtered.var_names]
        if existing_markers:
            sc.pl.umap(adata_filtered, color=existing_markers, 
                      save=f"_{sample_name}_{cell_type}_markers.pdf")
    
    # Find marker genes for clusters
    sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata_filtered, n_genes=25, sharey=False, 
                           save=f"_{sample_name}_markers.pdf")
    
    # Save processed data
    adata_filtered.write(os.path.join(OUTPUT_DIR, f"{sample_name}_processed.h5ad"))
    
    # Store for integration
    sample_adatas[sample_name] = adata_filtered
    all_samples.append(adata_filtered)

# Integration if we have multiple samples
if len(all_samples) > 1:
    print("\nIntegrating all samples...")
    
    # List of AnnData objects
    adata_list = list(sample_adatas.values())
    
    # Method 1: Basic concatenation for comparison
    adata_concat = adata_list[0].concatenate(
        adata_list[1:], batch_key='sample_name', index_unique='-'
    )
    
    # PCA and UMAP on concatenated data
    sc.pp.neighbors(adata_concat)
    sc.tl.umap(adata_concat)
    sc.tl.leiden(adata_concat, resolution=0.5)
    
    # Plot to check batch effects
    sc.pl.umap(adata_concat, color=['leiden', 'sample_name', 'tissue_type', 'tumor_stage'], 
              wspace=0.4, save="_concat_umap.pdf")
    
    # Method 2: Integration with Harmony (if installed)
    try:
        # Try to use Harmony for integration
        sc.external.pp.harmony_integrate(adata_concat, 'sample_name')
        sc.pp.neighbors(adata_concat, use_rep='X_pca_harmony')
        sc.tl.umap(adata_concat)
        sc.tl.leiden(adata_concat, resolution=0.5)
        
        # Plot integrated results
        sc.pl.umap(adata_concat, color=['leiden', 'sample_name', 'tissue_type', 'tumor_stage'], 
                  wspace=0.4, save="_harmony_umap.pdf")
        
        # Split by sample to check integration
        sc.pl.umap(adata_concat, color='leiden', split_by='sample_name', 
                  save="_harmony_split_umap.pdf")
        
        print("Harmony integration complete")
    except ImportError:
        print("Harmony not available, using basic concatenation only")
    
    # Method 3: Alternative integration using BBKNN (if that doesn't work)
    try:
        import bbknn
        # Reset the original concatenated object
        adata_bbknn = adata_list[0].concatenate(
            adata_list[1:], batch_key='sample_name', index_unique='-'
        )
        
        # Rerun PCA
        sc.tl.pca(adata_bbknn)
        
        # Run BBKNN
        bbknn.bbknn(adata_bbknn, batch_key='sample_name')
        sc.tl.umap(adata_bbknn)
        sc.tl.leiden(adata_bbknn, resolution=0.5)
        
        # Plot BBKNN integration
        sc.pl.umap(adata_bbknn, color=['leiden', 'sample_name', 'tissue_type', 'tumor_stage'], 
                  wspace=0.4, save="_bbknn_umap.pdf")
        
        # Save integrated object
        adata_bbknn.write(os.path.join(OUTPUT_DIR, "integrated_bbknn.h5ad"))
        print("BBKNN integration complete")
    except ImportError:
        print("BBKNN not available, skipping this integration method")
    
    # Save integrated dataset
    adata_concat.write(os.path.join(OUTPUT_DIR, "integrated_dataset.h5ad"))
    
    # Find marker genes for integrated clusters
    sc.tl.rank_genes_groups(adata_concat, 'leiden', method='wilcoxon')
    
    # Save markers to CSV
    markers_df = sc.get.rank_genes_groups_df(adata_concat, group=None)
    markers_df.to_csv(os.path.join(OUTPUT_DIR, "integrated_markers.csv"))
    
    print("Integration complete")

print(f"Analysis complete. Results saved to {OUTPUT_DIR}")