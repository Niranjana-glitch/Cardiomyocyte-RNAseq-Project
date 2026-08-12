# Scripts

## Overview

This folder contains the R scripts used for the complete RNA-seq analysis workflow of human induced pluripotent stem cell (hiPSC)-derived cardiomyocyte differentiation.

## Analysis Performed

The scripts perform:

* Import and preprocessing of RNA-seq count data
* Sample metadata integration
* Differential gene expression analysis using DESeq2
* Variance Stabilizing Transformation (VST)
* Principal Component Analysis (PCA)
* Sample correlation analysis
* Heatmap generation of top differentially expressed genes
* Gene Ontology (GO) enrichment analysis
* KEGG pathway enrichment analysis
* Volcano plot visualization
* Cardiac marker gene expression analysis

## Software and Packages

* R (v4.6.1)
* DESeq2
* EnhancedVolcano
* pheatmap
* clusterProfiler
* org.Hs.eg.db
* ggplot2
* enrichplot
* dplyr

## Main Script

* Cardiac_RNAseq_Analysis.R

This script reproduces the complete analysis pipeline and generates all result files and figures included in this repository.

