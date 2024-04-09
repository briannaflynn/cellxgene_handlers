import os
import pandas as pd
import anndata as ad
import sys
from cxg_logger import *
import datetime as dt
import numpy as np
import gc

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
            
            # Determine whether to use X or raw.X and get the gene names accordingly
            if row['raw_data_location'] == 'X':
                raw_counts_matrix = adata.X
                gene_names = adata.var_names
            elif row['raw_data_location'] == 'raw.X':
                raw_counts_matrix = adata.raw.X
                gene_names = adata.raw.var_names

            # Check type for matrix first, then convert the raw counts matrix to df and get observation annotations for each row
            if isinstance(raw_counts_matrix, np.ndarray):
                print("The raw counts matrix is already a NumPy ndarray.")
                raw_counts_df = pd.DataFrame(raw_counts_matrix, index=adata.obs_names, columns=gene_names)
            else:
                print("The raw counts matrix is not a NumPy ndarray, convert first.")
                raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(), index=adata.obs_names, columns=gene_names)
            
            print('Unfiltered shape', raw_counts_df.shape)

            if gene_filter_list != None:
                matching_gene_columns = [x for x in raw_counts_df.columns.to_list() if x in gene_filter_list]
                print(f'The following genes provided by the gene filter list are not present in this dataset: {[x for x in gene_filter_list if x not in matching_gene_columns]}')

                raw_counts_df = raw_counts_df[matching_gene_columns]
                print('Filtered shape', raw_counts_df.shape)
                
            observations_df = adata.obs.copy()

            # reset indexes for merge
            observations_df.reset_index(inplace=True)
            raw_counts_df.reset_index(inplace=True)
            # sometimes nice to print for debugging purposes but commenting for now to reduce logfile sizes
            #print('\n\nRaw counts dataframe:', raw_counts_df)
            #print('\n\nObservations dataframe:', observations_df)
            
            try:
                # check if the index columns match
                if raw_counts_df['index'].equals(observations_df['index']):
                # merge observations_df columns into raw_counts_df
                    raw_counts_df = pd.merge(raw_counts_df, observations_df, on='index', how='left')
                    print('\n\nMerged dataframe\n\n', f'Columns: {list(raw_counts_df.columns)}\n', raw_counts_df)

                else:
                    print('\n\nObservations and raw counts dataframe index columns don\'t match, exporting observations separately')
                    obs_path = f"{out_path}{ID}_observations.csv"
                    observations_df.to_csv(obs_path)
                csv_path = f"{out_path}{ID}_raw_counts.csv"
                raw_counts_df.to_csv(csv_path)
                print(f"\nExported raw counts to CSV: {csv_path}")

            except KeyError as error:
                print(f'Checking for compatible shape {file_path}\n{error}')
                if raw_counts_df.shape[0] == observations_df.shape[0]:
                    raw_counts_df = pd.concat([raw_counts_df, observations_df], axis=1)
                    csv_path = f"{out_path}{ID}_raw_counts.csv"
                    raw_counts_df.to_csv(csv_path)
                    print(f"\nExported raw counts to CSV: {csv_path}")
                pass
            
            # force collect every iteration
            gc.collect()

        except MemoryError as error:
            print('File not processed, memory issue - try exporting this dataset on another system', error)
            continue 
            
### Run commands
# should move to main.py with a better CLI
csv_request = sys.argv[1]
data_path = sys.argv[2]
export_path = sys.argv[3]

# gene filter list so files aren't so massive
# this is hard coded for now but should be updated for each user
ens_df = pd.read_csv('../ensembleid_FINAL.csv',header=None)
ens_list = ens_df[0].to_list()

df = pd.read_csv(csv_request)

print('#' * 25, 'STARTING EXPORT', '#' * 25)

date = str(dt.datetime.now())
date = date.replace(' ', '_')

logger = setup_logger('export', f'raw_data_export_{date}.log')

with StreamToLogger(logger, logging.INFO):
    try:
        export_raw_counts_to_csv(df, base_path_string=data_path, out_path=export_path, gene_filter_list=ens_list)
    except BaseException as error:
        print(f'BASE EXCEPTION: {error}, ID was NOT processed')
        pass