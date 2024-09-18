
# Quick Python scRNA-Seq Analysis Pipeline in Scanpy for Scoring Cell Types and Marker Genes of Interest

Overview

Do you wish to analyze scRNA-Seq data swiftly and efficiently, perhaps from the output of a pipeline such as nf-core/scRNA-seq?

This repository provides a Python script for analyzing single-cell RNA sequencing (scRNA-seq) data using Scanpy and related libraries. The pipeline encompasses data preprocessing, quality control, normalization, dimensionality reduction, clustering, and cell type annotation using scTYPE-inspired Python code.

Note: This script is not intended for absolute beginners. If you’re new to scRNA-Seq, consider taking a foundational course first, such as SIB Swiss Institute of Bioinformatics Training.

References

sc-type by Ianevski Aleksandr: https://github.com/IanevskiAleksandr/sc-type

Scanpy Documentation: https://scanpy.readthedocs.io/en/stable/installation.html

Requirements

Ensure you have a suitable environment to run the pipeline. I successfully ran it on an Mac M3 with 36 GB RAM, handling data from 20 patients and 120k cells at ease.

Required Files

1.	scRNA Data Files
	.loom files (e.g., PBMC_data_1.loom)
	For other formats (e.g., .h5ad), adapt the script accordingly.
2.	Annotations Files
	annotations.csv for each cell code with matching CellID if you intend to perform cell annotations.
3.	Gene List of Interest
	Example: human_marker_genes.csv for PBMC.

Contributing

Contributions are welcome! Please open an issue or submit a pull request for any enhancements or bug fixes.
