import os
import pandas as pd
import anndata as ad
import sys

def export_raw_counts_to_csv(df, base_path_string):
    """
    Exports raw count data from .h5ad files to CSV files based on a DataFrame input.
    
    Parameters:
    - df: DataFrame with columns 'collection_id', 'dataset_id', and 'raw_data_location'.
    - base_path_string: The base path to prepend to the paths constructed from the DataFrame.
    """
    # drop duplicate dataset ids so that no duplicate csv files get made
    print('\nStarting export\n\n')
    print(df)
    df_unique = df.drop_duplicates(subset='dataset_version_id', keep='first')
    print('\nDrop duplicate dataset IDs\n\n')
    print(df_unique)
    
    for idx, row in df.iterrows():

        # dataset version ID for naming, etc
        ID = row['dataset_version_id']

        print(ID)

        # construct the full path to the .h5ad file
        file_path = os.path.join(base_path_string, row['collection_id'], ID + '.h5ad')
        print('\nStarting with file:', file_path)

        if not os.path.exists(file_path):
            print(f'\n\nFile does not exist: {file_path}')
            continue

        adata = ad.read_h5ad(file_path)
        
        # Determine whether to use X or raw.X and get the gene names accordingly
        print('\n\nRaw data location:', row['raw_data_location'])
        if row['raw_data_location'] == 'X':
            raw_counts_matrix = adata.X
            gene_names = adata.var_names
            print(gene_names)
        elif row['raw_data_location'] == 'raw.X':
            raw_counts_matrix = adata.raw.X
            gene_names = adata.raw.var_names
            print(gene_names)
        
        # convert the raw counts matrix to df and get observation annotations for each row
        raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(), index=adata.obs_names, columns=gene_names)
        observations_df = adata.obs.copy()

        # reset indexes for merge
        observations_df.reset_index(inplace=True)
        raw_counts_df.reset_index(inplace=True)
        print('\n\nRaw counts dataframe:', raw_counts_df)
        print('\n\nObservations dataframe:', observations_df)
        
        # check if the index columns match
        if raw_counts_df['index'].equals(observations_df['index']):
            # merge observations_df columns into raw_counts_df
            merge_df = pd.merge(raw_counts_df, observations_df, on='index', how='left')
            print('\n\nMerged dataframe\n\n', f'Columns: {list(merge_df.columns)}\n', merge_df)
        else:
            print('\n\nObservations and raw counts dataframe index columns don\'t match, exporting observations separately')
            obs_path = f"{ID}_observations.csv"
            observations_df.to_csv(obs_path)
        
        csv_path = f"{ID}_raw_counts.csv"
        merge_df.to_csv(csv_path)
        print(f"\n\nExported raw counts to CSV: {csv_path}")

        sys.exit()

        

### example run

csv_request = sys.argv[1]
df = pd.read_csv(csv_request)
print('Starting test')
data = pd.read_csv('../data/aggregated_metadata_json_with_celltype_version.csv')[['dataset_id', 'dataset_version_id']].drop_duplicates(subset='dataset_version_id', keep='first')
print(data)
print(data.columns)
print(df)

merged_data = df.merge(data, on = 'dataset_id', how = 'inner')
print(merged_data)
print(merged_data.columns)
export_raw_counts_to_csv(merged_data, base_path_string="/work/projects/BioITeam/data/CELLxGENE/collections/")

###### the file is actually the dataset version ID
