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

def aggregate_and_modify_counts(df):
    """
    Aggregates numeric columns of a DataFrame by 'cell_type', computing mean and variance for each
    Uses pandas series 'value_counts()' method to get sample sizes for each cell type, and merges this into the aggregated df
    This is specifically useful for the combined variance calculation that will be performed after obtaining each cell type
    """
    
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    aggregations = {
        col: [('_mean', 'mean'), ('_var', 'var')]
        for col in numeric_cols
    }

    aggregated_data = df.groupby('cell_type', observed=False).agg(aggregations).reset_index()
    
    # Join MultiIndex column names into a single level, without this step, you won't be able to perform merges by column name 
    # '' allows you to combine without adding additional string
    aggregated_data.columns = [''.join(col).strip() for col in aggregated_data.columns.values]

    # merge sample size column of df.cell_type.value_counts dataframe into aggregated df
    aggregated_data = aggregated_data.merge(df.cell_type.value_counts().reset_index(name='sample_size'), on = 'cell_type')

    return aggregated_data

def generate_pandas_dataframe(file_path):
    """
    Generates a dataframe object using the raw count data from the anndata (.h5ad) object
    This data is accessible in the 'adata.raw.X', I presume that the data in 'adata.X' is this raw count data normalized in some way
    Specifically, this method filters the dataset down to only participants that are healthy (under 'disease' they are listed as 'normal') and adults
    """
    adata = ad.read_h5ad(file_path)

    observations = adata.obs.reset_index()
    print(observations.columns)
 
    # Check for columns containing 'age' and 'disease'
    age_columns = [col for col in observations.columns if 'age' in col.lower()]
    disease_columns = [col for col in observations.columns if 'disease' in col.lower()]

    # Ensure the assumption that columns contain 'age' and 'disease'
    if not age_columns:
        print("No columns containing 'age' were found.")
        return
    if not disease_columns:
        print("No columns containing 'disease' were found.")
        return

    # If columns are found, print them
    print(f"Columns containing 'age': {age_columns}")
    print(f"Columns containing 'disease': {disease_columns}")

    sys.exit()
    
    print(observations.development_stage.value_counts())
    # raw count data matrix
    raw_counts_matrix = adata.raw.X

    # convert the raw counts matrix to a pandas DataFrame - index will be the cell IDs (from adata.obs_names) and the columns will be the gene names (from adata.raw.var['feature_name'] or adata.raw.var_names)
    raw_counts_df = pd.DataFrame(raw_counts_matrix.toarray(),index=adata.obs_names, columns=adata.raw.var_names)
    raw_counts_df = raw_counts_df.reset_index()

    # subset the data to only healthy adults
    m_adult = observations[observations['age_group'] == 'adulthood']
    m_healthy_adult = m_adult[m_adult['disease'] == 'normal']
    healthy_adult_idx = m_healthy_adult['index'].to_list()
    filt_raw_counts_df = raw_counts_df.loc[raw_counts_df['index'].isin(healthy_adult_idx)]
    df = observations.merge(filt_raw_counts_df, on = 'index')

    # get cell type strings and replace spaces with underscores, making file and directory naming easier
    cell_types = [s.replace(" ", "_") for s in df.cell_type.unique()]

    aggregated_data = aggregate_and_modify_counts(df)
    fname = file_path.split('/')[-1]
    if not os.path.exists('../data/aggregated'):
        os.makedirs('../data/aggregated')

    aggregated_data.to_csv('../data/aggregated/aggregated_' + fname[:-5] + '.csv') 
    
    # just for visualization
    ei = ("#" * 25) + " "
    ej = " " + ("#" * 25)

    print(f'\n{ei}'+f'aggregated data for file: {file_path}'+f'{ej}\n')
    print(aggregated_data)

    for c in df.cell_type.unique():
        print(f'\n{ei}{c}{ej}\n')
        ag_c = aggregated_data[aggregated_data['cell_type'] == c]
        cell = c.replace(' ', '_')
        dir_name = '../data/' + cell 
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        print(ag_c)
        
        ag_c.to_csv(dir_name + '/' + fname + '_' + cell + '.csv')
        
def process_h5ad(file_path, summarize_file_structure=False):
    """
    Function for exploaratory analysis of anndata file structure
    """
    with h5py.File(file_path, 'r') as f:
        
        if summarize_file_structure:
            print("File structure of: ", file_path)
            f.visit(print)

        if 'raw/X/data' in f:
            print(file_path)
            print('\nRaw data available, generating tables\n')
            generate_pandas_dataframe(file_path)
        else:
            pass

def read_single_col_file(list_file):
    
    with open(list_file, 'r') as file:
        lines = [line.strip() for line in file.readlines()]

    return lines


