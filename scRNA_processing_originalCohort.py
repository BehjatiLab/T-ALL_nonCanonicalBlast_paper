########## T-ALL -- 061 Scanpy combined ########## 

##### Set up

# Load libraries
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc

# Print session info
sc.logging.print_header()

# Set working directory
#os.chdir('/lustre/scratch126/casm/team274sb/bl10/T-ALL/')



### Define variables

# Define s_genes and g2m_genes
s_genes =  ["MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "MRPL36", "E2F8"]
g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]



### Load metadata

# Load manifest
manifest = pd.read_csv("Data/TALL_manifest.csv")
manifest = manifest[manifest['Experiment'] == "GEX"]
manifest = manifest[manifest['Omit'] == "No"]
manifest.sort_values(by = ['Sample_ID'], inplace = True, ignore_index = True)

# Generate list_of_sample_ID
list_of_sample_ID = manifest['Sample_ID']
print(list_of_sample_ID)
print("Number of items =", len(list_of_sample_ID))

# Load qc_limits
qc_limits = pd.read_csv("Data/TALL_qc_limits.csv")





##### Load individual adata objects

list_of_adata = []

for i, sample_ID in enumerate(list_of_sample_ID):
    
    print(i)
    print(sample_ID)
    
    # Read SoupX adjusted count matrix
    adata = sc.read_10x_mtx("Intermediate/011_SoupX/"+sample_ID+"/soupX_matrix", var_names = 'gene_symbols', make_unique = True)
    adata.obs_names = [sample_ID+"::"+x for x in adata.obs_names]
    
    # Set sample_ID and cell_ID
    adata.obs['sample_ID'] = sample_ID
    adata.obs['cell_ID'] = adata.obs_names
    
    # Calculate QC metrics including n_genes, n_counts, percent_mito and percent_ribo
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo'], percent_top = None, log1p = False, inplace = True)
    adata.obs = adata.obs.rename(columns = {
        'n_genes_by_counts' : 'n_genes', 
        'total_counts' : 'n_counts', 
        'pct_counts_mito' : 'percent_mito', 
        'pct_counts_ribo' : 'percent_ribo'})
    
    # Input doublet_score and predicted_doublet
    scrublet_df = pd.read_csv("Intermediate/012_Scrublet/"+sample_ID+".csv")
    mapping_dict = dict(zip(scrublet_df['barcode'], scrublet_df['doublet_score']))
    adata.obs['doublet_score'] = adata.obs['cell_ID'].map(mapping_dict)
    mapping_dict = dict(zip(scrublet_df['barcode'], scrublet_df['predicted_doublet']))
    adata.obs['predicted_doublet'] = adata.obs['cell_ID'].map(mapping_dict)
    
    # Define QC limits
    n_genes_min = float(qc_limits[qc_limits['sample_ID'] == sample_ID]['nFeature_RNA_min'].iloc[0])
    n_counts_min = float(qc_limits[qc_limits['sample_ID'] == sample_ID]['nCount_RNA_min'].iloc[0])
    percent_mito_max = float(qc_limits[qc_limits['sample_ID'] == sample_ID]['percent_mt_max'].iloc[0])
    
    # Filter good quality cells
    adata = adata[adata.obs['n_genes'] > n_genes_min]
    adata = adata[adata.obs['n_counts'] > n_counts_min]
    adata = adata[adata.obs['percent_mito'] < percent_mito_max]
    adata = adata[adata.obs['predicted_doublet'] == False]

    # Store in list_of_adata
    list_of_adata.append(adata)
    del adata
    


    
    
##### Combine all adata objects and perform preprocessing

# Concatenate adata objects
adata = anndata.concat(list_of_adata)

# Save adata object containing raw counts
adata.write("Intermediate/061_Scanpy_combined/TALL_combined_raw_counts.h5ad")

# Log-normalize
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

# Select highly variable genes
sc.pp.highly_variable_genes(adata)

# Set .raw to the log-normalized counts for gene expression analysis (needed as scaling will modify adata.X)
adata.raw = adata

# Scale
sc.pp.scale(adata)

# Run PCA
sc.tl.pca(adata, n_comps = 100)

# Compute the neighborhood graph
sc.pp.neighbors(adata, n_pcs = 100, n_neighbors = 50)

# Embed the neighborhood graph
sc.tl.umap(adata, min_dist = 0.5)

# Cluster the neighborhood graph with varying resolutions
for i, res in enumerate([1.0, 2.0, 3.0, 4.0, 5.0]):
    sc.tl.leiden(adata, resolution = res)
    adata.obs['leiden_'+str(res)] = adata.obs['leiden']

# Call back log-normalized counts in .raw
adata = adata.raw.to_adata()

# Cell cycle phase
sc.tl.score_genes_cell_cycle(adata, s_genes = s_genes, g2m_genes = g2m_genes)
adata.obs.rename(columns = {'phase' : 'cell_cycle_phase'}, inplace = True)

# Save adata.obs and adata
adata.obs.to_csv("Intermediate/061_Scanpy_combined/TALL_combined_obs.csv", index = False)
adata.write("Intermediate/061_Scanpy_combined/TALL_combined.h5ad")




