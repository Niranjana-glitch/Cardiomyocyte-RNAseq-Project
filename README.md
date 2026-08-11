# Cardiac Differentiation RNA-Seq Analysis

## Project Overview

This project analyzes RNA-Seq data from human pluripotent stem cell-derived cardiomyocytes at different differentiation stages (Day 4, Day 5, and Day 8) to identify transcriptional changes associated with cardiac development.

## Objectives

- Perform quality-controlled transcriptomic analysis.
- Identify differentially expressed genes (DEGs) between differentiation stages.
- Explore global expression patterns using PCA and clustering.
- Investigate biological processes associated with cardiac differentiation through Gene Ontology (GO) enrichment analysis.
- Examine the expression of known cardiac marker genes.

## Dataset

- Source: NCBI Gene Expression Omnibus (GEO)
- Platform: RNA-Seq
- Samples: Day 4, Day 5, and Day 8 cardiac differentiation stages
- Organism: Homo sapiens

## Workflow

### 1. Data Acquisition
- Download raw count matrix and sample metadata.
- Import data into R.

### 2. Differential Expression Analysis
- DESeq2 normalization
- Dispersion estimation
- Wald statistical testing
- Multiple testing correction (Benjamini–Hochberg)

### 3. Exploratory Data Analysis
- Principal Component Analysis (PCA)
- Sample-to-sample distance heatmap
- Hierarchical clustering

### 4. Functional Enrichment Analysis
- Gene Ontology (GO) Biological Process enrichment
- Visualization using dot plots and bar plots

### 5. Marker Gene Analysis
- Expression profiling of cardiac marker genes
- Comparative expression across differentiation stages

## Software and Packages

### R Version
- R 4.6.1

### Packages
- DESeq2
- ggplot2
- pheatmap
- clusterProfiler
- enrichplot
- org.Hs.eg.db
- dplyr
- tidyr
- EnhancedVolcano
- factoextra

## Results

### Principal Component Analysis
PCA revealed distinct transcriptional profiles across differentiation stages, indicating temporal changes in gene expression during cardiomyocyte development.

### Hierarchical Clustering
Samples clustered according to differentiation stage, demonstrating biological consistency among replicates.

### Differential Gene Expression
Significant DEGs were identified between differentiation stages, highlighting genes involved in cardiac development and cellular differentiation.

### GO Enrichment Analysis
Enriched biological processes included:

- Cell fate commitment
- Cell fate specification
- Gland development
- Endocrine system development
- Tissue morphogenesis
- Pattern specification process
- Regulation of hormone levels

### Cardiac Marker Analysis
Known cardiac marker genes exhibited stage-specific expression patterns consistent with progressive cardiomyocyte differentiation.

## Figures

- PCA Plot
- Sample Distance Heatmap
- DEG Heatmap
- GO Enrichment Dot Plot
- GO Enrichment Bar Plot
- Cardiac Marker Expression Plot

## Output Files

- Significant_DEGs.csv
- PCA_plot.png
- SampleDistanceHeatmap.png
- DEG_Heatmap.png
- GO_Dotplot.png
- GO_Barplot.png
- CardiacMarkerExpression.png

## Author

Shameem Fathima

B.Tech Biotechnology

Vel Tech High Tech Dr. Rangarajan Dr. Sakunthala Engineering College

Research Interests:
- Bioinformatics
- Transcriptomics
- Cardiovascular Biology
- Functional Genomics