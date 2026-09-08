# Multi-Omics Analysis Pipeline: PTEN and PTK2B Signaling in TNBC

## Overview
This repository contains the complete computational biology pipeline developed for my Master of Research (MRes) dissertation. The master script integrates multi-omics data to investigate longitudinal adaptations to PTEN reintroduction and PTK2B (Pyk2) signaling dynamics in Triple-Negative Breast Cancer (TNBC).

The pipeline processes raw bulk transcriptomics and mass spectrometry proteomics, reconstructs causal signaling networks, and validates these mechanistic findings using patient-derived single-cell RNA sequencing and clinical cohorts.

## Pipeline Architecture
The workflow is consolidated into a single master R script (`PTK2B_Master_Pipeline.R`) designed for sequential, linear execution. Major analytical modules include:

* **Pre-Processing & Imputation:** Bulk RNA-seq processing via `DESeq2` and LFQ proteomics analysis with missing-value imputation via `proDA`.
* **Dimensionality Reduction & Clustering:** PCA and MOFA+ multi-omics integration, followed by kinetic trajectory clustering using `Mfuzz`.
* **Activity Inference:** Transcription factor and pathway footprinting using `CollecTRI`, `RcisTarget`, and `PROGENy`.
* **Causal Network Reconstruction:** PKN filtering via `OmnipathR` and mechanistic network optimization using the `CARNIVAL` pipeline (Gurobi solver).
* **Clinical & Immune Correlations:** Integration of TCGA basal cohorts via `GSVA` and immune deconvolution mapping using `CIBERSORTx`.
* **Single-Cell Validation:** Spatial footprinting and doublet removal in TNBC patient atlas datasets (GSE176078) using `Seurat` and `harmony`.

## Prerequisites
To execute the pipeline, ensure the following core R packages are installed:
* **Core Data Sci:** `dplyr`, `tidyr`, `ggplot2`, `ComplexHeatmap`, `patchwork`
* **Omics Processing:** `DESeq2`, `proDA`, `Seurat`, `limma`, `GSVA`
* **Network & Systems Bio:** `CARNIVAL`, `decoupleR`, `progeny`, `dorothea`, `OmnipathR`, `RcisTarget`
* **Optimization:** Requires a valid local installation and license for the Gurobi solver (used for CARNIVAL).

## Data Availability
Due to file size constraints, raw `.mtx` count matrices, TCGA clinical `.txt` files, and raw mass spectrometry outputs are not hosted in this repository. Ensure that all raw data inputs and external metadata files are structured in their respective sub-directories (e.g., `/RNAseq_1/`, `/proteomics/`, `/TCGA+cibersort/`) relative to the working directory before executing the script.

## Author
* **Ben [Your Last Name]** 
* University of Southampton
