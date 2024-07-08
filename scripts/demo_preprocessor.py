from typing import Dict, Optional, Union, Tuple, List
import numpy as np
from scipy.sparse import issparse
import scanpy as sc
from scanpy.get import _get_obs_rep, _set_obs_rep
from anndata import AnnData
from datetime import datetime
import os
import sys

# add the parent directory "scripts" to the Python path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

#from cxg_logger import *
from adata_preprocessor import Preprocessor
from synthetic_data_generator import *

# example execution scripts for the Preprocessor class

# synthetic_adata = generate_synthetic_anndata(num_cells=1000, num_genes=1000, sparsity = 0.7, observations={'cell_type':['neuron'] * 1000, 'disease':['normal']* 1000, 'organism': ['human'] * 700 + ['mouse'] * 300})
# print('Synthetic data before observation filter')
# print(synthetic_adata.obs)

# filter_obs_preprocessor = Preprocessor(
#     use_key=None,
#     filter_gene_by_counts=False,
#     gene_filter=False,
#     filter_cell_by_counts=False,
#     normalize_total=False,
#     result_normed_key="X_normed",
#     log1p=False,
#     result_log1p_key="X_log1p",
#     subset_hvg=False,
#     hvg_use_key=None,
#     hvg_flavor="seurat_v3",
#     binning=None,
#     result_binned_key="X_binned",
#     even_binning=False,
#     filter_observations={'organism':'human'},
#     execute_filter_genes='off',
#     execute_filter_cells='off',
#     execute_normalize_total='off',
#     execute_log1p='off',
#     execute_subset_hvg='off',
#     execute_binning='off'
# )

# filter_obs_preprocessor(synthetic_adata)

# print('Synthetic data after observation filter via preprocessor class')
# # see 'filter_observations' and provided dict
# print(synthetic_adata.obs)

# """
# # Output
# Synthetic data before observation filter
#       cell_id cell_type disease organism
# 0      cell_0    neuron  normal    human
# 1      cell_1    neuron  normal    human
# 2      cell_2    neuron  normal    human
# 3      cell_3    neuron  normal    human
# 4      cell_4    neuron  normal    human
# ..        ...       ...     ...      ...
# 995  cell_995    neuron  normal    mouse
# 996  cell_996    neuron  normal    mouse
# 997  cell_997    neuron  normal    mouse
# 998  cell_998    neuron  normal    mouse
# 999  cell_999    neuron  normal    mouse

# [1000 rows x 4 columns]
# Synthetic data after observation filter via preprocessor class
#       cell_id cell_type disease organism
# 0      cell_0    neuron  normal    human
# 1      cell_1    neuron  normal    human
# 2      cell_2    neuron  normal    human
# 3      cell_3    neuron  normal    human
# 4      cell_4    neuron  normal    human
# ..        ...       ...     ...      ...
# 695  cell_695    neuron  normal    human
# 696  cell_696    neuron  normal    human
# 697  cell_697    neuron  normal    human
# 698  cell_698    neuron  normal    human
# 699  cell_699    neuron  normal    human

# [700 rows x 4 columns]
# """

# gene_filter_preprocessor = Preprocessor(
#     use_key=None,
#     filter_gene_by_counts=False,
#     gene_filter={"gene_id" : ["gene_0", "gene_1", "gene_2"]},
#     filter_cell_by_counts=False,
#     normalize_total=False,
#     result_normed_key="X_normed",
#     log1p=False,
#     result_log1p_key="X_log1p",
#     subset_hvg=False,
#     hvg_use_key=None,
#     hvg_flavor="seurat_v3",
#     binning=None,
#     result_binned_key="X_binned",
#     even_binning=False,
#     execute_filter_genes='on',
#     execute_filter_cells='off',
#     execute_normalize_total='off',
#     execute_log1p='off',
#     execute_subset_hvg='off',
#     execute_binning='off'
# )

# gene_filter_preprocessor(synthetic_adata)

# print(synthetic_adata.var)
# """
#   gene_id
# 0  gene_0
# 1  gene_1
# 2  gene_2
# """

# # lite version of class initialization
# # genes to keep assumes you're providing a list of indexes that correspond to 'var_names', print adata.var_names to check if index is gene names
# # print(synthetic_adata.var_names)

# lite_gene_filter_preprocessor = Preprocessor(
#     use_key=None,
#     genes_to_keep=["0", "1"],
#     execute_filter_genes='on'
# )

# lite_gene_filter_preprocessor(synthetic_adata)

# print(synthetic_adata.var)

# """
#   gene_id
# 0  gene_0
# 1  gene_1
# """

# synthetic_adata = generate_synthetic_anndata(num_cells=1000, num_genes=1000, sparsity = 0.7, observations={'cell_type':['neuron'] * 1000, 'disease':['normal']* 1000, 'organism': ['human'] * 700 + ['mouse'] * 300})

# raw_lite_gene_filter_preprocessor = Preprocessor(
#     use_key=None,
#     convert_raw=True,
#     genes_to_keep=["1"],
#     execute_filter_genes='on'
# )

# raw_lite_gene_filter_preprocessor(synthetic_adata)

# print(synthetic_adata.raw.var)

# print(synthetic_adata.var)

# """
# 2024-07-08 08:57:55,561 - INFO - Converting raw data matrix anndata.raw.X to main matrix anndata.X
# 2024-07-08 08:57:55,562 - INFO - Updating raw.var and raw.var_names to var and var_names
# 2024-07-08 08:57:55,562 - INFO - Filtering genes to keep only the specified genes ... 
#       gene_id
# 0      gene_0
# 1      gene_1
# 2      gene_2
# 3      gene_3
# 4      gene_4
# ..        ...
# 995  gene_995
# 996  gene_996
# 997  gene_997
# 998  gene_998
# 999  gene_999

# [1000 rows x 1 columns]
#   gene_id
# 1  gene_1
# """
