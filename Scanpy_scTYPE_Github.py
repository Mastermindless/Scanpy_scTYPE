# Simple install in your jupyter notebook
# !pip install scanpy 
# !pip install openpyxl
# !pip install loompy
# !pip install igraph
# !pip install leidenalg
# !pip install scipy

# scRNAseq_analysis.py

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import leidenalg
import scipy

# Set the working directory to the script's location for relative path handling
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# List of Loom files in dir files:
loom_files = [
    "files/PBMC_data_1.loom",
    "files/PBMC_data_2.loom",
	"files/PBMC_data_final.loom"
]

# Initialize an empty list to store AnnData objects
adatas = []

# Process each Loom file
for file in loom_files:
    adata = sc.read_loom(file)

    # Check and ensure unique gene names
    duplicates_before = adata.var_names.duplicated().sum()
    if duplicates_before > 0:
        print(f"[{file}] Duplicated gene names before making unique: {duplicates_before}")
    else:
        print(f"[{file}] No duplicated gene names before making unique.")

    adata.var_names_make_unique()

    duplicates_after = adata.var_names.duplicated().sum()
    if duplicates_after > 0:
        print(f"[{file}] Duplicated gene names after making unique: {duplicates_after}")
    else:
        print(f"[{file}] All gene names are unique after applying .var_names_make_unique()")

    # Add batch information using the file name
    adata.obs['batch'] = os.path.basename(file)
    adatas.append(adata)

# Concatenate all AnnData objects into one
combined_adata = sc.concat(adatas, join='outer', label='batch', keys=loom_files)

# Clean up
del adatas, adata

# Export genes and barcodes
combined_adata.var_names.to_series().to_csv("output/genes.csv")
combined_adata.obs_names.to_series().to_csv("output/barcodes.csv")

# Annotate single cells
annotations = pd.read_csv("data/annotations.csv", index_col=0)

# Join annotations with the combined AnnData object
combined_adata.obs = combined_adata.obs.join(annotations)

# Export specific genes to check for available markers of interest
genes_of_interest = ['PTPRC', 'CD8A']
for gene in genes_of_interest:
    gene_level = combined_adata[:, gene].X.toarray().flatten()
    pd.DataFrame({f'{gene}_lvl': gene_level}, index=combined_adata.obs_names).to_csv(f"output/{gene}_lvl_data.csv")
    combined_adata.obs[f'{gene}_lvl'] = gene_level

# Quality control - calculate mitochondrial gene percentage
combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Normalize and log-transform the data
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(combined_adata, flavor='seurat', n_top_genes=2000)

# Scale the data
sc.pp.scale(combined_adata, max_value=10)

# Perform PCA
sc.tl.pca(combined_adata, svd_solver='arpack')

# Visualize PCA
sc.pl.pca(combined_adata, color='sample', legend_loc='right margin', save='_patient.png')

# Elbow plot to determine the number of PCs
sc.pl.pca_variance_ratio(combined_adata, log=True, n_pcs=50, save='_elbow_plot.png')

# Compute the neighborhood graph
sc.pp.neighbors(combined_adata, n_pcs=42)

# Perform Leiden clustering
sc.tl.leiden(
    combined_adata,
    resolution=0.8,
    flavor='igraph',       # Specify the backend as 'igraph'
    n_iterations=2,       # Number of iterations
    directed=False        # Ensure the graph is undirected
)

# Compute UMAP for visualization
sc.tl.umap(combined_adata)
sc.pl.umap(combined_adata, color=['leiden', 'sample'], legend_loc='on data', save='_check_map_clusters.png')

# High-resolution exports
visualization_colors = ['leiden', 'sample', 'STAGE', 'TMB', 'Gender', 'Response']
for color in visualization_colors:
    sc.pl.umap(combined_adata, color=[color], legend_loc='right margin', save=f'_{color}_umap_clusters.pdf')

# Implement scTYPE logic
# Read marker genes from an Excel file
marker_db = pd.read_excel("data/marker_genes.xlsx")
tissue_type = 'scTYPE'  # Adjust as necessary
cell_markers = marker_db[marker_db['tissueType'] == tissue_type].copy()

# Function to process gene symbols
def process_genes(gene_series):
    gene_series = gene_series.str.replace(' ', '').str.replace('///', ',')
    gene_lists = gene_series.apply(lambda x: [gene.upper() for gene in x.split(',') if gene and gene.upper() != 'NA'] if pd.notnull(x) else [])
    return gene_lists

cell_markers['gs_positive'] = process_genes(cell_markers['geneSymbolmore1'])
cell_markers['gs_negative'] = process_genes(cell_markers['geneSymbolmore2'])

# Create dictionaries for positive and negative gene sets
gs_positive = dict(zip(cell_markers['cellName'], cell_markers['gs_positive']))
gs_negative = dict(zip(cell_markers['cellName'], cell_markers['gs_negative']))

# Prepare the scaled data for scoring
if scipy.sparse.issparse(combined_adata.X):
    scRNAseqData = pd.DataFrame(
        combined_adata.X.todense().T,
        index=combined_adata.var_names.str.upper(),
        columns=combined_adata.obs_names
    )
else:
    scRNAseqData = pd.DataFrame(
        combined_adata.X.T,
        index=combined_adata.var_names.str.upper(),
        columns=combined_adata.obs_names
    )

# Function to compute scTYPE scores
def sctype_score(scRNAseqData, gs_positive, gs_negative):
    scores = pd.DataFrame(index=gs_positive.keys(), columns=scRNAseqData.columns)

    for cell_type, pos_genes in gs_positive.items():
        neg_genes = gs_negative.get(cell_type, [])

        pos_genes_present = [gene for gene in pos_genes if gene in scRNAseqData.index]
        neg_genes_present = [gene for gene in neg_genes if gene in scRNAseqData.index]

        pos_score = scRNAseqData.loc[pos_genes_present].mean(axis=0) if pos_genes_present else 0
        neg_score = scRNAseqData.loc[neg_genes_present].mean(axis=0) if neg_genes_present else 0

        scores.loc[cell_type] = pos_score - neg_score

    return scores

# Compute scTYPE scores
scores = sctype_score(scRNAseqData, gs_positive, gs_negative)

# Assign cell types per cluster
clusters = combined_adata.obs['leiden']
unique_clusters = clusters.unique()
cluster_results = []

for cl in unique_clusters:
    cells_in_cl = clusters[clusters == cl].index
    es_max_cl = scores[cells_in_cl].sum(axis=1).sort_values(ascending=False)
    total_cells = len(cells_in_cl)
    top_cell_type = es_max_cl.idxmax()
    top_score = es_max_cl.max()
    cluster_results.append({
        'cluster': cl,
        'type': top_cell_type,
        'score': top_score,
        'ncells': total_cells
    })

cluster_results_df = pd.DataFrame(cluster_results)

# Set low-confidence clusters to 'Unknown'
cluster_results_df.loc[
    cluster_results_df['score'] < cluster_results_df['ncells'] / 4, 'type'
] = 'Unknown'

# Map cell types to the AnnData object
combined_adata.obs['scTYPE'] = combined_adata.obs['leiden'].map(
    cluster_results_df.set_index('cluster')['type']
)

# Export cluster annotations
cluster_results_df.to_csv("output/cluster_annotations.csv", index=False)
combined_adata.obs.to_csv("output/single_cell_annotations.csv")

# High-resolution UMAP export for scTYPE and gene levels
sc.pl.umap(combined_adata, color=['scTYPE'], legend_loc='right margin', save='_scTYPE_umap_clusters.pdf')
for gene in genes_of_interest:
    sc.pl.umap(combined_adata, color=[f'{gene}_lvl'], legend_loc='right margin', color_map='viridis_r', save=f'_{gene}_umap_clusters.pdf')


