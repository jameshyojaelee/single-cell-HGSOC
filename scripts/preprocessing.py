#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Single-cell RNA-seq preprocessing script for HGSOC data using Scanpy.
This script processes CellRanger outputs and merges technical replicates 
based on the GSM sample names from the metadata.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import defaultdict

# Set random seed for reproducibility
np.random.seed(42)

# Set plotting style
sc.settings.set_figure_params(dpi=100, frameon=False)
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
plt.rcParams['figure.figsize'] = (8, 8)

# Paths
METADATA_PATH = "/gpfs/commons/home/jameslee/HGSOC/metadata/GSE184880_metadata.csv"
CELLRANGER_DIR = "/gpfs/commons/home/jameslee/HGSOC/output/cellranger"
OUTPUT_DIR = "/gpfs/commons/home/jameslee/HGSOC/output/scanpy"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load metadata
metadata = pd.read_csv(METADATA_PATH)
print(f"Loaded metadata with {len(metadata)} entries")

# Group SRR runs by GSM sample name
sample_groups = defaultdict(list)
for _, row in metadata.iterrows():
    sample_groups[row['Sample Name']].append(row['Run'])

print(f"Found {len(sample_groups)} unique GSM samples")

# Process each GSM sample
all_samples = []

for gsm_name, srr_list in sample_groups.items():
    print(f"\nProcessing {gsm_name} with {len(srr_list)} sequencing runs")
    
    # Check if all SRR data exists
    srr_adatas = []
    for srr_id in srr_list:
        h5_path = os.path.join(CELLRANGER_DIR, srr_id, "outs", "filtered_feature_bc_matrix.h5")
        if not os.path.exists(h5_path):
            print(f"Warning: CellRanger output for {srr_id} not found at {h5_path}")
            continue
        
        # Load data
        print(f"Loading {srr_id}...")
        adata = sc.read_10x_h5(h5_path)
        
        # Add metadata
        meta_row = metadata[metadata['Run'] == srr_id].iloc[0]
        adata.obs['run_id'] = srr_id
        adata.obs['sample_name'] = gsm_name
        adata.obs['age'] = meta_row['age']
        adata.obs['pathology'] = meta_row['pathology']
        adata.obs['tissue_type'] = meta_row['tissue_type']
        adata.obs['tumor_stage'] = meta_row['tumor_stage']
        
        # Add a unique barcode suffix to avoid conflicts when merging
        adata.obs_names = [f"{bc}_{srr_id}" for bc in adata.obs_names]
        
        # Basic QC metrics
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        
        srr_adatas.append(adata)
    
    if not srr_adatas:
        print(f"No data found for {gsm_name}, skipping")
        continue
    
    # Merge all sequencing runs for this sample
    if len(srr_adatas) > 1:
        combined = srr_adatas[0].concatenate(srr_adatas[1:], batch_key='sequencing_run')
        print(f"Merged {len(srr_adatas)} runs for {gsm_name}, resulting in {combined.n_obs} cells")
    else:
        combined = srr_adatas[0]
        combined.obs['sequencing_run'] = '0'
        print(f"Single run for {gsm_name}, {combined.n_obs} cells")
    
    # Basic preprocessing for this sample
    print(f"Preprocessing {gsm_name}...")
    
    # Filter cells based on QC metrics (adjust thresholds as needed)
    sc.pp.filter_cells(combined, min_genes=200)
    sc.pp.filter_genes(combined, min_cells=3)
    
    # Filter out cells with high mitochondrial gene fraction
    combined.var['mt'] = combined.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], inplace=True)
    combined = combined[combined.obs['pct_counts_mt'] < 20].copy()
    
    print(f"After QC: {combined.n_obs} cells, {combined.n_vars} genes")
    
    # Save the processed individual sample
    combined.write(os.path.join(OUTPUT_DIR, f"{gsm_name}_processed.h5ad"))
    
    # Store for integration later
    all_samples.append(combined)

# Integrate all samples
if len(all_samples) > 0:
    print("\nIntegrating all samples...")
    integrated = all_samples[0].concatenate(all_samples[1:], batch_key='sample_name')
    print(f"Integrated dataset: {integrated.n_obs} cells, {integrated.n_vars} genes")
    
    # Normalize and log transform
    sc.pp.normalize_total(integrated, target_sum=1e4)
    sc.pp.log1p(integrated)
    
    # Find variable genes
    sc.pp.highly_variable_genes(integrated, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='sample_name')
    print(f"Found {sum(integrated.var.highly_variable)} highly variable genes")
    
    # Scale data
    sc.pp.scale(integrated, max_value=10)
    
    # Run PCA
    sc.tl.pca(integrated, svd_solver='arpack')
    
    # Batch correction using harmony (if installed) or bbknn
    try:
        import harmonypy
        sc.external.pp.harmony_integrate(integrated, 'sample_name')
        use_rep = 'X_pca_harmony'
    except ImportError:
        print("Harmony not available, using standard PCA")
        use_rep = 'X_pca'
    
    # Compute neighborhood graph
    sc.pp.neighbors(integrated, n_neighbors=10, n_pcs=30, use_rep=use_rep)
    
    # Run UMAP and clustering
    sc.tl.umap(integrated)
    sc.tl.leiden(integrated, resolution=0.5)
    
    # Calculate marker genes
    sc.tl.rank_genes_groups(integrated, 'leiden', method='wilcoxon')
    
    # Plot results
    sc.pl.umap(integrated, color=['leiden', 'sample_name', 'tissue_type', 'tumor_stage'], 
               ncols=2, save=f"_umap_clusters_samples.png")
    
    # Create separate UMAPs for specific groupings
    sc.pl.umap(integrated, color='tissue_type', save="_umap_tissue_type.png")
    sc.pl.umap(integrated, color='tumor_stage', save="_umap_tumor_stage.png")
    
    # Save the integrated dataset
    integrated.write(os.path.join(OUTPUT_DIR, "integrated_dataset.h5ad"))
    
    # Save marker genes to CSV
    result = integrated.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    marker_df = pd.DataFrame({
        group + '_' + key: result[key][group]
        for group in groups
        for key in ['names', 'pvals', 'pvals_adj', 'scores']
    })
    
    marker_df.to_csv(os.path.join(OUTPUT_DIR, "marker_genes.csv"))
    
    print(f"Analysis complete. Results saved to {OUTPUT_DIR}")
else:
    print("No samples were processed successfully!")