# Quick python scRNA-Seq analysis pipeline in Scanpy and scTYPE

You want to analyze scRNASeq data from e.g. the output of a pipline such as https://nf-co.re/scrnaseq/2.6.0/

Check out references:
https://github.com/IanevskiAleksandr/sc-type
https://scanpy.readthedocs.io/en/stable/installation.html

This repository contains a Python script for analyzing single-cell RNA sequencing (scRNA-seq) data using Scanpy and related libraries. The pipeline includes data preprocessing, quality control, normalization, dimensionality reduction, clustering, and cell type annotation using scTYPE inspired python code. This is not a pure beginners script so maybe take a course first if you are new to scRNASeq: https://www.sib.swiss/training/course/2021-03-sc-Transcriptomics

You need:

1. scRNA_datat.loom files 

2. annotations.csv files for each cell code with machin CellID if you want to do cell annotations

3. Gene list of interests. Here is an example list for human_marker_genes for PBMC.

Skript is ready to test for PBMC_data_1.loom data For other formats (e.g .hd5 etc.), please adapt script accordingly. 

LGM
