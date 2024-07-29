import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
from collections import defaultdict
from copy import copy
import re

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

list_of_sample_ID = ["P058_D0", "P058_D28"]

# Load adata for P058
adata = sc.read_h5ad("Data/P058.h5ad")



### Figure 4E: ZBTB16 expression on UMAP of P058 blasts (same as Figure 1E)

fig, ax = plt.subplots(figsize = (4.4, 3.5))
reds_cmap = copy(matplotlib.cm.Reds)
reds_cmap.set_under("lightgrey")
sc.pl.umap(adata, color = 'ZBTB16', color_map = reds_cmap, vmin = 0.00001, vmax = 4.4, ax = ax)
plt.savefig("Plots/Fig4E_P058_UMAP_ZBTB16.pdf", bbox_inches = 'tight')
plt.close()



### Figure 4F: TCR gene usage on UMAP of P058 blasts (called by TRUST4)

# Concatenate multiple vdj objects from the same patient
list_of_vdj = []
for j, sample_ID in enumerate(list_of_sample_ID): 
    vdj = pd.read_csv("Data/"+sample_ID+"_barcode_report.tsv", sep = "\t")
    vdj = vdj.rename(columns = {'#barcode' : 'barcode'})
    vdj['barcode'] = sample_ID + '::' + vdj['barcode']
    list_of_vdj.append(vdj)
vdj = pd.concat(list_of_vdj)
del list_of_vdj
    
# Annotate V/D/J/C calls and nucleotide sequence
vdj[['chain1_V', 'chain1_D', 'chain1_J', 'chain1_C', 'chain1_CDR3_seq']] = vdj['chain1'].str.split(pat = ',', expand = True).iloc[:, 0:5]
vdj[['chain2_V', 'chain2_D', 'chain2_J', 'chain2_C', 'chain2_CDR3_seq']] = vdj['chain2'].str.split(pat = ',', expand = True).iloc[:, 0:5]
    
# Annotate clones by V/J calls and nucleotide sequence
vdj['chain1_VJ_clone'] = vdj['chain1_V'] + ',' + vdj['chain1_J']
vdj['chain1_VJseq_clone'] = vdj['chain1_V'] + ',' + vdj['chain1_J'] + ',' + vdj['chain1_CDR3_seq']
vdj['chain2_VJ_clone'] = vdj['chain2_V'] + ',' + vdj['chain2_J']
vdj['chain2_VJseq_clone'] = vdj['chain2_V'] + ',' + vdj['chain2_J'] + ',' + vdj['chain2_CDR3_seq']

# Merge vdj with adata.obs
vdj = vdj.set_index('barcode')
adata.obs = adata.obs.merge(vdj, how = 'left', left_index = True, right_index = True)

list_of_tcr_gene = ['TRBV4-2*01', 'TRBV4-1*01', 'TRBJ2-1*01', 'TRBJ1-2*01']
list_of_tcr_locus = ['chain1_V', 'chain1_V', 'chain1_J', 'chain1_J']
list_of_colour = ['#1D71B8', '#D62728', '#1D71B8', '#D62728']

for j, tcr_gene in enumerate(list_of_tcr_gene): 
    
    tcr_locus = list_of_tcr_locus[j]
    colour = list_of_colour[j]
    
    fig, ax = plt.subplots(figsize = (4, 3.5))
    colour_dict = defaultdict(lambda: '#D3D3D3')
    colour_dict[tcr_gene] = colour
    sc.pl.umap(adata, color = tcr_locus, groups = tcr_gene, palette = colour_dict, na_color = '#D3D3D3', na_in_legend = False, size = 20, ax = ax)
    plt.savefig("Plots/Fig4F_P058_TCR_TRUST4_UMAP_"+tcr_gene.replace('*', '_')+".pdf", bbox_inches = 'tight')
    plt.close()

    

### Figure 4G: Copy-Number Alterations on UMAP of P058 blasts

# Concatenate posterior probability for each sample
list_of_pp_all = []
for i, sample_ID in enumerate(list_of_sample_ID):
    pp_all = pd.read_csv("Data/"+sample_ID+"_scRNA_pp.csv")
    pp_all['cellID'] = pp_all['cellID'].str.replace(sample_ID+'_', sample_ID+'::')
    list_of_pp_all.append(pp_all)
pp_all = pd.concat(list_of_pp_all)
del list_of_pp_all

list_of_region = ['chr9p', 'chr17q']
list_of_colour = ['#000000', '#1D71B8']

for i, region in enumerate(list_of_region):
    
    colour = list_of_colour[i]
    
    # Subset posterior probability for only region of interest
    pp = pp_all[pp_all['regionID'] == region]

    # Transfer CN state and probability of CN state
    mapping_dict = dict(zip(pp['cellID'], pp['mostLikelyState']))
    adata.obs['CN_state'] = adata.obs['cell_ID'].map(mapping_dict)
    mapping_dict = dict(zip(pp['cellID'], pp['maxPostProb']))
    adata.obs['CN_state_prob'] = adata.obs['cell_ID'].map(mapping_dict)

    # Recode CN state as 'Aberrant', 'Normal' or 'Unsure' (probability < 0.95)
    adata.obs.loc[adata.obs['CN_state'] == 'abbFrac', 'CN_state'] = 'Aberrant'
    adata.obs.loc[adata.obs['CN_state'] == 'normFrac', 'CN_state'] = 'Normal'
    adata.obs.loc[adata.obs['CN_state_prob'] < 0.95, 'CN_state'] = 'Unsure'
    adata.obs.loc[adata.obs['CN_state'].isna(), 'CN_state'] = 'Unsure'
    
    fig, ax = plt.subplots(figsize = (4, 3.5))
    colour_dict = defaultdict(lambda: '#D3D3D3')
    colour_dict['Aberrant'] = colour
    sc.pl.umap(adata, color = 'CN_state', groups = ['Aberrant'], palette = colour_dict, title = region, na_color = '#D3D3D3', na_in_legend = False, size = 20, ax = ax)
    plt.savefig("Plots/Fig4G_P058_CN_altered_UMAP_"+region+".pdf", bbox_inches = 'tight')
    plt.close()



### Figure 4H: SNVs on UMAP of P058 blasts

list_of_cluster = ['Shared', 'CloneA']
list_of_colour = ['#000000', '#1D71B8']

for i, cluster in enumerate(list_of_cluster):
    
    colour = list_of_colour[i]
    
    # Load SNVs
    df = pd.read_csv("Data/P058_scRNA_SNV.csv")
    
    # Subset for mutation
    df = df[df['Cluster'] == cluster]
    
    # Transfer to adata
    adata.obs['Locus'] = 'Normal'
    adata.obs.loc[adata.obs['cell_ID'].isin(df['CB']), 'Locus'] = 'Mutated'
    
    # Skip if no cell in adata object contains the mutation
    if all(adata.obs['Locus'] == 'Normal'):
        continue

    # Plot UMAP (vector)
    fig, ax = plt.subplots(figsize = (4, 3.5))
    colour_dict = defaultdict(lambda: '#D3D3D3')
    colour_dict['Mutated'] = colour
    sc.pl.umap(adata, color = 'Locus', groups = ['Mutated'], palette = colour_dict, title = cluster, na_color = '#D3D3D3', na_in_legend = False, size = 30, ax = ax)
    plt.savefig("Plots/Fig4H_P058_SNV_UMAP_"+cluster+".pdf", bbox_inches = 'tight')
    plt.close()