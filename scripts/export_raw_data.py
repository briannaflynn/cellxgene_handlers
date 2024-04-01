import os
import pandas as pd
import anndata as ad

def export_raw_counts_to_csv(df, base_path_string):
    """
    Exports raw count data from .h5ad files to CSV files based on a DataFrame input.
    
    Parameters:
    - df: DataFrame with columns 'collection_id', 'dataset_id', and 'raw_data_location'.
    - base_path_string: The base path to prepend to the paths constructed from the DataFrame.
    """
    # drop duplicate dataset ids so that no duplicate csv files get made
    df_unique = df.drop_duplicates(subset='dataset_id', keep='first')
    
    for idx, row in df.iterrows():

        # dataset ID for naming, etc
        ID = row['dataset_id']
        # construct the full path to the .h5ad file
        file_path = os.path.join(base_path_string, row['collection_id'], ID + '.h5ad')
        
        if not os.path.exists(file_path):
            print(f"File does not exist: {file_path}")
            continue
        
        adata = ad.read_h5ad(file_path)
        
        # Determine whether to use X or raw.X and get the gene names accordingly
        if row['raw_data_location'] == 'X':
            raw_counts_matrix = adata.X
            gene_names = adata.var_names
        elif row['raw_data_location'] == 'raw.X':
            raw_counts_matrix = adata.raw.X
            gene_names = adata.raw.var_names
        
        # convert the raw counts matrix to df and get observation annotations for each row
        raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(), index=adata.obs_names, columns=gene_names)
        raw_counts_df = raw_counts_df.reset_index()
        observations_df = adata.obs.copy()

        # reset indexes for merge
        observations_df.reset_index(inplace=True)
        raw_counts_df.reset_index(inplace=True)
        
        # check if the index columns match
        if raw_counts_df['index'].equals(observations_df['index']):
            # merge observations_df columns into raw_counts_df
            raw_counts_df = pd.merge(raw_counts_df, observations_df, on='index', how='left')
        else:
            print('Observations and raw counts dataframe index columns don\'t match, exporting observations separately')
            obs_path = f"{ID}_observations.csv"
            observations_df.to_csv(obs_path)
        
        csv_path = f"{ID}_raw_counts.csv"
        raw_counts_df.to_csv(csv_path)
        print(f"Exported raw counts to CSV: {csv_path}")
