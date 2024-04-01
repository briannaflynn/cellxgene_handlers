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
    for idx, row in df.iterrows():
        # Construct the full path to the .h5ad file
        file_path = os.path.join(base_path_string, row['collection_id'], row['dataset_id'] + '.h5ad')
        
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
        
        # Convert the raw counts matrix to a pandas DataFrame
        raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(), index=adata.obs_names, columns=gene_names)
        raw_counts_df = raw_counts_df.reset_index()
        
        csv_path = f"{row['dataset_id']}_raw_counts.csv"
        raw_counts_df.to_csv(csv_path, index=False)
        print(f"Exported raw counts to CSV: {csv_path}")
