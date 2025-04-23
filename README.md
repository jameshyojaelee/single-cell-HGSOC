# Single-Cell Transcriptomic Analysis of High-Grade Serous Ovarian Cancer (HGSOC)

## Overview
This repository contains the analysis pipeline and results for studying transcriptomic profiles of high-grade serous ovarian cancer (HGSOC) at single-cell resolution. The project focuses on understanding immune evasion mechanisms and therapeutic resistance in HGSOC using scRNA-seq data.

## Research Questions
- How do transcriptomic profiles differ between healthy ovarian tissue and HGSOC cancer cells?
- How do HGSOC tumors change at the transcriptomic level in response to various immune checkpoint inhibitor (ICI) treatments?
- What are the genetic and phenotypic states of key immune cell populations within the tumor microenvironment (TME), and how do they contribute to immune evasion?

## Hypothesis
Distinct transcriptomic signatures and immune cell phenotypes within the TME of HGSOC contribute significantly to chemotherapy and ICI resistance. Specifically, differences in gene regulatory networks and intercellular communications among immune cells, cancer cells, and stromal components are correlated with varying treatment responses in patients.

## Data Sources

| GEO Accession | Study | Description |
|---------------|-------|-------------|
| GSE184880 | Xu, Fang, et al. | scRNA-seq of normal ovarian tissue vs. HGSOC ([Publication](https://aacrjournals.org/clincancerres/article/28/16/3590/707396/Single-Cell-RNA-Sequencing-Reveals-the-Tissue)) |
| GSE160755 | Wan, et al. | • scRNA-seq of untreated organoid and parent tumor<br>• Bulk RNA-seq of HGSOC Organoids with different combinations of ICI (PD-1/PD-L1) ([Publication](https://aacrjournals.org/cancerres/article/81/1/158/649480/Enhanced-Efficacy-of-Simultaneous-PD-1-and-PD-L1)) |

## Methods

### Preprocessing and Quality Control
- Alignment of raw fastq files using CellRanger v9.0.0 with GRCh38 human reference genome (including intronic reads)
- Quality control using Scanpy (Python) and Seurat (R)
- Doublet identification with Scrublet
- Batch correction using ComBat for integration across samples

### Analysis Pipeline
- Dimensionality reduction (PCA, UMAP) and clustering
- Differential gene expression analysis comparing:
  - Normal ovarian tissue vs. HGSOC samples
  - Untreated vs. ICI-treated tumors
- Immune cell annotation using ceLLamma and scType
- Cell-cell communication network analysis with CellChat/NicheNet
- Gene regulatory network analysis using scPRINT

## Expected Results
- Identification of differentially expressed genes and noncoding RNAs between normal and cancerous ovarian tissues
- Characterization of immune cell populations associated with therapeutic resistance
- Discovery of novel biological pathways and regulatory networks altered in HGSOC
- Identification of potential biomarkers or therapeutic targets

## Challenges and Solutions
- **Challenge**: Batch effects across different datasets  
  **Solution**: Applying robust computational methods for batch correction (Harmony, ComBat)

- **Challenge**: Accurate immune cell subtype annotation  
  **Solution**: Leveraging established marker gene sets, literature-based manual correction

- **Challenge**: Complex TME interpretation  
  **Solution**: Potential integration with spatial transcriptomics data (if available) 