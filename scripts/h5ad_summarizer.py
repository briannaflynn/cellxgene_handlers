import h5py
import sys
import anndata as ad
import pandas as pd
import scipy
import numpy as np
import os

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

def aggregate_and_modify_counts(df):
    """
    Aggregates numeric columns of a DataFrame by 'cell_type', computing mean, variance, and count for each.
    Keeps the first '_count' column, renames it to 'sample_size', and drops the other '_count' columns.
    """
    
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    aggregations = {
        col: [('_mean', 'mean'), ('_var', 'var')]
        for col in numeric_cols
    }

    aggregated_data = df.groupby('cell_type').agg(aggregations).reset_index()
    # Join MultiIndex column names into a single level
    aggregated_data.columns = [''.join(col).strip() for col in aggregated_data.columns.values]

    # get sample size from 
    aggregated_data = aggregated_data.merge(df.cell_type.value_counts().reset_index(name='sample_size'), on = 'cell_type')

    return aggregated_data


def generate_pandas_dataframe(file_path):

    adata = ad.read_h5ad(file_path)
    # observations only
    observations = adata.obs.reset_index()
    print(observations.columns)
    # Access the raw count data
    raw_counts_matrix = adata.raw.X

    # Convert the raw counts matrix to a pandas DataFrame
    # The index will be the cell IDs (from adata.obs_names) and the columns will be the gene names (from adata.raw.var['feature_name'] or adata.raw.var_names)
    raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(),index=adata.obs_names, columns=adata.raw.var_names)
    raw_counts_df = raw_counts_df.reset_index()

    # only healthy adults
    m_adult = observations[observations['age_group'] == 'adulthood']
    m_healthy_adult = m_adult[m_adult['disease'] == 'normal']
    healthy_adult_idx = m_healthy_adult['index'].to_list()

    filt_raw_counts_df = raw_counts_df.loc[raw_counts_df['index'].isin(healthy_adult_idx)]
    print(filt_raw_counts_df)
    
    merged_df = observations.merge(filt_raw_counts_df, on = 'index')

    print('make aggregated data')

    df = merged_df

    # get cell type strings and replace spaces with underscores
    cell_types = [s.replace(" ", "_") for s in df.cell_type.unique()]

    aggregated_data = aggregate_and_modify_counts(df)

    fname = file_path.split('/')[-1]
    aggregated_data.to_csv('../data/aggregated_' + fname[:-5] + '.csv') 

    for c in df.cell_type.unique():
        print(f'\n{c}\n')
        ag_c = aggregated_data[aggregated_data['cell_type'] == c]
 
        cell = c.replace(' ', '_')
        dir_name = '../data/' + cell 
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        print(ag_c)
        
        ag_c.to_csv(dir_name + '/' + fname + '_' + cell + '.csv')
        

if __name__ == '__main__':
  file_path = sys.argv[1] # .hda5 file is the first argument
  #summarize_h5ad(file_path)
  generate_pandas_dataframe(file_path)


###### TODO #########
# write bash program to merge all text files into one, make sure that headers are the same 
# compute mean of means