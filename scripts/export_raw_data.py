import os
import pandas as pd
import anndata as ad
import sys
from cxg_logger import *
import datetime as dt
import numpy as np
import gc
from adata_preprocessor import Preprocessor

def export_raw_counts_to_csv(df, base_path_string:str, out_path:str, exclude_IDs:list = None, gene_filter_list=None):
    """
    Exports raw count data from .h5ad files to CSV files based on a DataFrame input.
    
    Parameters:
    - df: DataFrame with columns 'collection_id', 'dataset_id', and 'raw_data_location'.
    - base_path_string: The base path to prepend to the paths constructed from the DataFrame.
    """
    # drop duplicate dataset ids so that no duplicate csv files get made
    print('\nStarting export\n\n')
    df_unique = df.drop_duplicates(subset='dataset_version_id', keep='first')

    # on or off for running gene filter from Preprocessor
    if gene_filter_list != None:
        run_gene_filter = 'on'
    else:
        run_gene_filter = 'off'

    # if theres a list of dataset version IDs to exclude, remove these first 
    if exclude_IDs != None:
        exclude = exclude_IDs
        df_unique = df_unique[~df_unique['dataset_version_id'].isin(exclude)]

    for idx, row in df_unique.iterrows():
        try:
            # dataset version ID for naming, etc
            ID = row['dataset_version_id']

            print(f"Dataset version ID being processed: {ID}")

            # construct the full path to the .h5ad file
            file_path = os.path.join(base_path_string, row['collection_id'], ID + '.h5ad')
            print('Starting with file:', file_path)

            if not os.path.exists(file_path):
                print(f'File does not exist: {file_path}')
                continue

            # read h5ad file as anndata object
            adata = ad.read_h5ad(file_path)
            print("Original gene list length:", len(adata.var_names))

            # Determine whether to use X or raw.X and get the gene names accordingly
            if row['raw_data_location'] == 'X':
                gene_filter_preprocessor = Preprocessor(
                    # False if X
                    convert_raw = False,
                    genes_to_keep = gene_filter_list,
                    execute_filter_genes=run_gene_filter
                )
                gene_filter_preprocessor(adata)
                
            elif row['raw_data_location'] == 'raw.X':
                gene_filter_preprocessor = Preprocessor(
                    # True if raw.X
                    convert_raw = True,
                    genes_to_keep = gene_filter_list,
                    execute_filter_genes=run_gene_filter
                )
                gene_filter_preprocessor(adata)

            # if raw, convert_raw True converts adata.raw.X to adata.X
            raw_counts_matrix = adata.X
            gene_names = adata.var_names
            print("Filtered gene list length:", len(gene_names))

            # Check type for matrix first, then convert the raw counts matrix to df and get observation annotations for each row
            if isinstance(raw_counts_matrix, np.ndarray):
                print("The raw counts matrix is already a NumPy ndarray.")
                raw_counts_df = pd.DataFrame(raw_counts_matrix, index=adata.obs_names, columns=gene_names)
            else:
                print("The raw counts matrix is not a NumPy ndarray, convert first.")
                raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(), index=adata.obs_names, columns=gene_names)
                
            observations_df = adata.obs.copy()    
            obs_path = f"{out_path}{ID}_observations.csv"
            observations_df.to_csv(obs_path)
            csv_path = f"{out_path}{ID}_raw_counts.csv"
            raw_counts_df.to_csv(csv_path)
            if raw_counts_df.shape[0] == observations_df.shape[0]:
                print("Observation and raw counts shapes are compatible")
            print(f"\nExported raw counts and observations to CSV: {csv_path}")
            # force collect every iteration
            gc.collect()

        except MemoryError as error:
            print('File not processed, memory issue - try exporting this dataset on another system', error)
            continue 
