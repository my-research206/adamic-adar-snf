# Adamic-Adar Similarity Network Fusion

This repository contains the R implementation of the Adamic-Adar refinement of Similarity Network Fusion (AA-SNF) used for multi-omics patient similarity analysis.

## Contents

- `adamic_adar_snf_analysis_breast_cancer_cells.R` — R implementation of the AA-SNF analysis and evaluation pipeline.

## Dataset

The analysis uses the GSE22210 breast cancer dataset obtained from the NCBI Gene Expression Omnibus (GEO).

## Method

AA-SNF refines the initial SNF patient-affinity network using Adamic-Adar local neighborhood topology before the standard SNF diffusion and fusion process.

## Software

The implementation uses selected functions from the SNFtool R package.
