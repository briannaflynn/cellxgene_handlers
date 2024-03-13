import h5py
import sys
import anndata as ad
import pandas as pd
import scipy

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
    #raw_counts_df = raw_counts_df.reset_index()
    raw_counts_df = raw_counts_df.astype({col: pd.SparseDtype("float", 0) for col in raw_counts_df.columns})
    print(raw_counts_df.reset_index())
    #merged_df = observations[['index', 'cell_type', 'age_group', 'cell_type_ontology_term_id']].merge(raw_counts_df, on = 'index')
    #print(merged_df)
    #print(merged_df.columns)


if __name__ == '__main__':
  file_path = sys.argv[1] # .hda5 file is the first argument
  #summarize_h5ad(file_path)
  generate_pandas_dataframe(file_path)
