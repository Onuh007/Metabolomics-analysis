# Metabolomics Analysis Pipeline

## Overview
A complete R-based metabolomics analysis pipeline comparing  metabolite profiles between Control and Disease groups.

## Methods
- Quality control and missing value imputation
- Log2 transformation and autoscaling
- Univariate statistical testing (t-test, Wilcoxon)
- Multiple testing correction (Benjamini-Hochberg FDR)
- Fold change and log2FC calculation
- Principal Component Analysis (PCA)
- PLS-DA with VIP scoring
- Volcano plots, box plots, correlation heatmaps

## Requirements
- R version 4.5
- Packages: tidyverse, FactoMineR, factoextra, mixOmics, corrplot

## Key Findings
- All 10 metabolites showed statistically significant elevation in the Disease group after BH correction (p < 0.05).
- PC1 explained 93.9% of total variance with complete group separation observed.

## Author
Augustine Chukwu Onuh Ph.D
Specialisation: Metabolomics, Bioinformatics  
Languages: R, Python

## Analysis Report
Full interactive report available here:
https://rpubs.com/Onuh007/1438987
