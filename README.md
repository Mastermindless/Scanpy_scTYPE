# Quick python scRNA-Seq analysis pipeline in Scanpy and scTYPE

You want to analyze scRNASeq data from e.g. the output of a pipline such as https://nf-co.re/scrnaseq/2.6.0/

Check out references:
https://github.com/IanevskiAleksandr/sc-type
https://scanpy.readthedocs.io/en/stable/installation.html

This repository contains a Python script for analyzing single-cell RNA sequencing (scRNA-seq) data using Scanpy and related libraries. The pipeline includes data preprocessing, quality control, normalization, dimensionality reduction, clustering, and cell type annotation using scTYPE inspired python code. This is not a pure beginners script so maybe take a course first if you are new to scRNASeq: https://www.sib.swiss/training/course/2021-03-sc-Transcriptomics

You need:

A decent notebook. I ran it on a M3 with 36Gb with 20 patients and 120k cells.

1. scRNA_datat.loom files 

2. annotations.csv files for each cell code with machin CellID if you want to do cell annotations

3. Gene list of interests. Here is an example list for human_marker_genes for PBMC.

Skript is ready to test for PBMC_data_1.loom data For other formats (e.g .hd5 etc.), please adapt script accordingly. 

LGM


scanpy==1.9.1
pandas==1.5.3
numpy==1.24.4
matplotlib==3.7.2
seaborn==0.12.2
leidenalg==0.8.4
scipy==1.10.1
loompy==3.0.7
openpyxl==3.1.2
igraph==0.10.7


