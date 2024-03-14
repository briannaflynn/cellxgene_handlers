import h5py
import sys
import anndata as ad
import pandas as pd
import scipy
import numpy as np

'''Usage
path to your h5ad file is supplied as argument
i.e.
python h5ad_summarizer.py my_dataset.h5ad
'''

def summarize_h5ad(file_path):
    with h5py.File(file_path, 'r') as f:
        print("File structure of:", file_path)
        f.visit(print)  # Print the hierarchy of the file

        # If the file follows the AnnData structure, it will have a 'X' dataset under some group that represents the main data matrix
        print("\nSummary of 'X' dataset:")
        data = f['X']['data']
        #print(data)
        print(f"Shape: {data.shape}")
        print(f"Data type: {data.dtype}")

        # layers
        print('layers')
        print(f['layers'])

        # basic statistics of dense data
        if isinstance(data, h5py.Dataset):
            # this loads entire dataset, be careful if low on memory
            data_loaded = data[:] 
            print(data_loaded)
            print(f"Mean: {data_loaded.mean()}")
            print(f"Standard deviation: {data_loaded.std()}")
                

        # print out the names of variables (genes, features) if they exist
        if 'var' in f:
            print('\nFirst 10 variable names in feature_name/categories')
            print(f['var/feature_name/categories'][:])
            print('\nFirst 10 gene IDs')
            print(f['var/gene_ids'][:10])
            genes = f['var']['gene_ids']
            # Shape of gene_ids (19492,)
            print(f'Shape of gene_ids {genes.shape}')

        if 'obs' in f:
            print('Cell type')
            print(f['obs/cell_type/categories'].shape)
            print('n count RNA example')
            print(f['obs/nCount_RNA'][:10])
            print(f['obs/cell_type_ontology_term_id/codes'][:])
            print(f['obs/cell_type/codes'][:])
            print(f['obs']['cell_type']['codes'].shape)

def check_if_sparse(file_path):

    adata = ad.read_h5ad(file_path)

    print('Get adata layers')
    print(adata.layers.keys())
    is_sparse= False
    if isinstance(adata.X, (scipy.sparse.csv_matrix, scipy.sparse.csc_matrix)):
        is_sparse = True
    print(f"Sparse matrix detected? {is_sparse}")

import pandas as pd

def calculate_combined_std(df):
    """
    Calculate the combined standard deviation from multiple groups.
    
    Parameters:
    - df: pandas DataFrame with columns 'mean', 'variance', and 'sample_size'.
    
    Returns:
    - combined_std: Combined standard deviation of the groups.
    """
    # Calculate the weighted average of the means (grand mean)
    grand_mean = (df['mean'] * df['sample_size']).sum() / df['sample_size'].sum()

    # Calculate weighted sum of variances and the correction for the mean differences
    weighted_sum_of_variances = ((df['variance'] * (df['sample_size'] - 1)).sum() +
                                 (df['sample_size'] * ((df['mean'] - grand_mean) ** 2)).sum())

    # Calculate the combined variance
    combined_variance = weighted_sum_of_variances / (df['sample_size'].sum() - len(df))

    # Calculate the combined standard deviation
    combined_std = combined_variance ** 0.5

    return combined_std

def generate_pandas_dataframe(file_path):

    adata = ad.read_h5ad(file_path)
    # print('adata layers')
    # print(adata.layers.keys())
    # if scipy.sparse.issparse(adata.X):
    #     data_matrix = adata.X.toarray()
    # else:
    #     data_matrix = adata.X
    # print(adata.obs.columns)
    # df = pd.DataFrame(data_matrix, columns=adata.var_names, index=adata.obs_names)

    # print(df)

    # observations only
    observations = adata.obs.reset_index()
    print(observations)
    # Access the raw count data
    raw_counts_matrix = adata.raw.X

    # Convert the raw counts matrix to a pandas DataFrame
    # The index will be the cell IDs (from adata.obs_names) and the columns will be the gene names (from adata.raw.var['feature_name'] or adata.raw.var_names)
    raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(),index=adata.obs_names, columns=adata.raw.var_names)
    
    raw_counts_df = raw_counts_df.astype({col: pd.SparseDtype("float", 0) for col in raw_counts_df.columns})
    raw_counts_df = raw_counts_df.reset_index()
    
    merged_df = observations.merge(raw_counts_df, on = 'index')
    print('merge dataframe')
    print(merged_df)
    print(merged_df.columns)

    # subset to just healthy adults 
    m_adult = merged_df[merged_df['age_group'] == 'adulthood']
    m_healthy_adult = m_adult[m_adult['disease'] == 'normal']
    print('make aggregated data')
    aggregated_data = m_healthy_adult.groupby('cell_type').agg({col: ['group_mean', 'variance', 'sample_size'] for col in m_healthy_adult.select_dtypes(include=[np.number]).columns})
    print(aggregated_data)
    aggregated_data.to_csv(file_path[:-5] + '.csv')    

if __name__ == '__main__':
  file_path = sys.argv[1] # .hda5 file is the first argument
  #summarize_h5ad(file_path)
  generate_pandas_dataframe(file_path)
