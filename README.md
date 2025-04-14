# MeERCAT: Metabolite and Expressed RNA Cross-modality Analysis Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) <!-- Choose your license and update badge -->

A Python package designed to process, integrate, and analyze paired transcriptomics (RNA-Seq) and metabolomics (mass spectrometry) data. MeERCAT facilitates a workflow from raw data loading and cleaning to correlation analysis (Spearman) and dimensionality reduction using Non-negative Matrix Factorization (NMF) for identifying joint patterns.

## Overview

Mass spectrometry data integration presents numerous challenges that make it difficult to incorporate into many sequencing-based multi-omics workflows. This package provides tools to handle common challenges in integrating RNA-Seq count data and metabolomics abundance data, specifically focusing on:

*   Loading data from multiple experimental batches/files.
*   Cleaning and standardizing feature names and sample identifiers.
*   Creating robust composite sample IDs based on experimental metadata.
*   Aligning datasets to common samples.
*   Filtering features based on variance.
*   Performing pairwise correlation analysis (Spearman) with p-value adjustment (Benjamini-Hochberg FDR).
*   Applying concatenated Non-negative Matrix Factorization to identify joint latent factors across both omics types.
*   Evaluating NMF results using standard metrics.
*   Generating visualizations like Volcano plots, heatmaps, and scatter plots.

## Features

*   Loads metadata and multiple RNA-Seq/Metabolite CSV files.
*   Extracts experimental IDs and metadata from filenames or data columns.
*   Cleans feature names (e.g., removing gene versions).
*   Handles duplicate features by summing.
*   Cleans metabolite data (numeric conversion, NaN handling).
*   Filters out QC/Blank samples based on sample names.
*   Creates consistent composite sample identifiers.
*   Aligns RNA and Metabolite datasets to common samples.
*   Filters low-count RNA genes.
*   Normalizes RNA data using CPM + log2(x+1).
*   Filters features by variance (Top N).
*   Calculates Spearman correlations between all filtered gene-metabolite pairs.
*   Adjusts p-values using Benjamini-Hochberg FDR.
*   Performs concatenated NMF after scaling data (MaxAbsScaler) and imputing NaNs (Zero Imputation).
*   Calculates NMF evaluation metrics (Reconstruction Error, Component Contribution, Sparsity, H-Correlation).
*   Generates various plots: Rho distribution, Volcano plot, Global clustermap, Top feature heatmap, Top correlation scatter plots, NMF H-correlation heatmap, NMF k-evaluation plot.
*   Saves intermediate and final results to structured output directories.
*   Configurable parameters via `meercat/config.py`.

## Installation

Currently, this package is intended for local installation from source.

1.  **Clone the repository:**
    ```bash
    git clone <your-repo-url>
    cd meercat
    ```
2.  **Create and activate a virtual environment (Recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```
3.  **Install dependencies and the package:**
    ```bash
    pip install --upgrade pip build
    pip install .
    ```
    This command uses the `pyproject.toml` file to install the package and its dependencies listed there.

## Input Data Requirements

The package expects input data to be organized in a specific structure within a base input directory. By default, it looks for an `input_data` subdirectory relative to where the analysis script is run.