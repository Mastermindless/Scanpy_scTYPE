# scRNA-seq Analysis Pipeline

Check out:
https://github.com/IanevskiAleksandr/sc-type
https://scanpy.readthedocs.io/en/stable/installation.html

This repository contains a Python script for analyzing single-cell RNA sequencing (scRNA-seq) data using Scanpy and related libraries. The pipeline includes data preprocessing, quality control, normalization, dimensionality reduction, clustering, and cell type annotation using scTYPE.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Data](#data)
- [Output](#output)
- [Contributing](#contributing)
- [License](#license)

## Installation

1. **Clone the repository:**

    ```bash
    git clone https://github.com/yourusername/scRNAseq_analysis.git
    cd scRNAseq_analysis
    ```

2. **Set up a virtual environment (optional but recommended):**

    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3. **Install the required packages:**

    ```bash
    pip install -r requirements.txt
    ```

## Usage

1. **Prepare your data:**

    - Place all Loom files in the `files/` directory.
    - Ensure that annotation files (`annotations.csv` and `marker_genes.xlsx`) are in the `data/` directory.

2. **Run the analysis script:**

    ```bash
    python scRNAseq_analysis.py
    ```

3. **Review Outputs:**

    - All output files will be saved in the `output/` directory.
    - Visualizations will be saved alongside the script or in specified subdirectories.

## Directory Structure
