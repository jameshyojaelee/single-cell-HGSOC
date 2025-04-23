#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set working directory
os.chdir("/gpfs/commons/home/jameslee/HGSOC/output")

# Read the data
count_matrix = pd.read_csv("count_matrix.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)
umap_coords = pd.read_csv("umap_coordinates.csv", index_col=0)

# Create AnnData object
adata = sc.AnnData(X=count_matrix.T)  # Transpose to get cells as rows
adata.obs = metadata
adata.obsm['X_umap'] = umap_coords.values

# Basic preprocessing
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_qc_metrics.pdf')

# Normalize and log transform
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='_highly_variable_genes.pdf')

# Scale the data
# sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_pca_variance.pdf')

# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP (using existing coordinates)
sc.pl.umap(adata, color=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
           save='_umap_qc.pdf')

# Clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color='leiden', save='_umap_clusters.pdf')

# Find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes.pdf')

# Save the processed data
adata.write('processed_data.h5ad')

# Print some basic statistics
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Number of clusters: {len(adata.obs['leiden'].unique())}") 

# Assume cell types are annotated in adata.obs['cell_type']
# Example annotation (adjust markers as needed)
adata.obs['cell_type'] = 'Other'  # Default
adata.obs.loc[adata[:, ['TCF7', 'ITGAE']].X.mean(axis=1) > 1, 'cell_type'] = 'TRMstem CD8'
adata.obs.loc[adata[:, ['PDCD1', 'LAG3']].X.mean(axis=1) > 1, 'cell_type'] = 'Exhausted CD8'
adata.obs.loc[adata[:, ['CD68']].X > 1, 'cell_type'] = 'Macrophage'

# Plot UMAP with all cell types
plt.figure(figsize=(8, 6))
sc.pl.umap(adata, color='cell_type', title='UMAP of Immune Cell Types in Ovarian Cancer TME', show=False)
plt.savefig('umap_all_cell_types.png')
plt.close()

# Highlight specific cell types
plt.figure(figsize=(8, 6))
sc.pl.umap(adata, color='cell_type', groups=['TRMstem CD8', 'Exhausted CD8'], 
           title='UMAP Highlighting TRMstem and Exhausted CD8+ Cells', show=False)
plt.savefig('umap_highlighted_cell_types.png')
plt.close()