# sc-RNA-Seq

# Single Cell RNA-seq Analysis of PBMC Data

This repository contains an R-based analysis pipeline for single cell RNA sequencing (scRNA-seq) data. The pipeline processes and analyzes a dataset of 2,700 single cells (sequenced on the Illumina NextSeq 500) using popular bioinformatics tools and packages such as Seurat and SingleR. The analysis covers quality control, normalization, variable feature identification, dimensionality reduction (PCA and UMAP), clustering, differential expression, and cell type annotation.

## Overview

The provided Jupyter Notebook (using the R kernel) implements the following tasks:

- **Loading Required Libraries:** Utilizes packages such as `dplyr`, `Seurat`, `patchwork`, `ggplot2`, `SingleR`, and others for a seamless analysis workflow.
- **Data Import:** Reads in the filtered gene expression matrix generated from 10X Genomics data.
- **Quality Control & Filtering:** Calculates mitochondrial gene percentages and creates QC plots (violin and scatter plots). Filters cells based on minimum feature counts and mitochondrial content.
- **Normalization & Feature Selection:** Normalizes the data and identifies highly variable features.
- **Dimensionality Reduction & Clustering:** Scales the data, runs PCA for dimension reduction, performs clustering (using graph-based methods), and visualizes clusters with UMAP.
- **Differential Expression Analysis:** Identifies differentially expressed markers across clusters and extracts the top 3 markers for each cluster.
- **Cell Type Annotation:** Uses the SingleR package with the Human Primary Cell Atlas reference to assign cell type labels and compares them to the Seurat-derived clusters.

## Repo Contents
. ├── filtered_gene_bc_matrices/ │ └── hg19/ # Directory containing raw 10X Genomics data. ├── sc_analysis_results_pbmc/ # Output directory for generated plots and results. ├── analysis_notebook.ipynb # Jupyter Notebook with the analysis pipeline. └── README.md # This file.

## Analysis Workflow

The analysis pipeline covers the following steps:

### Data Loading and Library Setup
- **Required libraries are loaded.**
- **Data Import:** Data from the 10X Genomics filtered matrix is imported and a Seurat object is created.

### Quality Control (QC)
- **Mitochondrial Gene Calculation:** Mitochondrial gene percentages are calculated.
- **Plotting:** Violin and scatter plots (for `nFeature_RNA`, `nCount_RNA`, and `percent.mt`) are generated.
- **Cell Filtering:** Cells are filtered based on defined quality metrics.

### Normalization & Variable Feature Identification
- **Normalization:** The data is normalized using the LogNormalize method.
- **Feature Identification:** Highly variable genes are identified using the `FindVariableFeatures` function.

### Dimensionality Reduction
- **Scaling & PCA:** Data scaling and PCA computation are performed.
- **Diagnostic Plots:** Diagnostic plots such as PCA feature plots and elbow plots are generated.

### Clustering and UMAP Visualization
- **Clustering:** Clusters are identified via graph-based algorithms (`FindNeighbors` and `FindClusters`).
- **UMAP Visualization:** UMAP is used to visualize the cell clusters.

### Differential Expression Analysis
- **Marker Identification:** Differentially expressed marker genes are identified using the `FindAllMarkers` function.
- **Marker Extraction:** Top marker genes per cluster are extracted and saved.

### Cell Type Annotation with SingleR
- **Annotation:** The SingleR package is used with the Human Primary Cell Atlas reference dataset to annotate cell types.
- **Comparison & Visualization:** UMAP visualizations, score heatmaps, and contingency tables comparing SingleR results to Seurat clusters are generated.

## Data Description

- **Input:** The dataset consists of 2,700 single cells from PBMC samples, sequenced using the Illumina NextSeq 500. Data is provided as a 10X Genomics filtered gene expression matrix.

- **Output:** The analysis outputs include:
  - Quality control plots (e.g., `QC_violin.png`, `QC_scatter.png`)
  - Variable feature plots (e.g., `VariableFeatures.png`)
  - PCA diagnostics (e.g., `PCA_plot.png`, `Elbow_plot.png`, `PCA_top_features.txt`)
  - UMAP visualizations (e.g., `UMAP_clusters.png`, `UMAP_clusters_labeled.png`, `UMAP_SingleR_labels.png`)
  - Differential expression results (`top3_markers_per_cluster.csv`)
  - Contingency tables and heatmaps (`cluster_vs_SingleR_table.csv`, `cluster_vs_SingleR_heatmap.png`)

## Prerequisites

- **Jupyter Notebook with an R Kernel:** Ensure that you have Jupyter Notebook installed and configured to use R.
- **R (version 4.x or higher recommended)**
- **R Packages:**
  - CRAN packages: `dplyr`, `patchwork`, `ggplot2`, `pheatmap`, `DT`
  - Bioconductor packages: `Seurat`, `SingleR`, `celldex`

