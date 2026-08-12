# Raw Count Data

## Dataset Description

This folder contains the raw RNA-seq count matrix and associated metadata used in the analysis.

## Files

### GSE116574_count_matrix.csv

Raw gene expression count matrix downloaded from the Gene Expression Omnibus (GEO).

### metadata.csv

Experimental metadata describing sample groups and conditions.

## Study Design

Samples represent different stages of cardiomyocyte differentiation:

* Day 4 (3D4)
* Day 5 (3D5)
* Day 8 (3D8)

## Data Source

GEO Accession:
GSE116574

Database:
NCBI Gene Expression Omnibus (GEO)

## Usage

These files serve as the primary input for:

* Differential expression analysis
* PCA
* Clustering
* Heatmap generation
* Functional enrichment analysis
* Marker gene analysis

The raw count matrix can be directly imported into DESeq2 for reproducible RNA-seq analysis.

