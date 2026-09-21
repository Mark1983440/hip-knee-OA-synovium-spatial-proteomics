# Hip–Knee OA Spatial Proteomics
Analysis scripts associated with the study:

**Spatial Proteomics Identifies Distinct Synovial Microenvironments in Hip and Knee Osteoarthritis**

## Overview

GeoMx DSP spatial proteomics data were analysed in R using workflows adapted from the `standR` package.

## Analysis workflow

The analysis includes:

- quality control
- assessment of candidate normalization factors
- TMM normalization and RUVg correction
- relative log expression (RLE) and principal component analysis (PCA) assessment
- differential protein abundance analysis using limma-voom
- hip-versus-knee comparisons within:
  - Intima CD45+
  - Intima CD45−
  - Subintima CD45+
  - Subintima CD45−
- generation of volcano plots and IPA input files


## Software

Analyses were performed in R version 4.5.1.

Main packages include:

- standR 1.12.0
- limma 3.64.1
- edgeR 4.6.3
- SpatialExperiment
- scater

## References

1. Tan CW, Berrell N, Donovan ML, Monkman J, Lawler C, Sadeghirad H, et al. The development of a high-plex spatial proteomic methodology for the characterisation of the head and neck tumour microenvironment. *NPJ Precis Oncol.* 2025;9(1):191.

2. Liu N, Bhuva DD, Mohamed A, Bokelund M, Kulasinghe A, Tan CW, Davis MJ. standR: spatial transcriptomic analysis for GeoMx DSP data. *Nucleic Acids Res.* 2024;52(1):e2.
