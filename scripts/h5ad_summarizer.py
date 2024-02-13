import h5py
import sys

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
            print(f"Shape: {data.shape}")
            print(f"Data type: {data.dtype}")

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
            print(f['var/feature_name/categories'][:10])
            print('\nFirst 10 gene IDs')
            print(f['var/gene_id'][:10])

if __name__ == '__main__':
  file_path = sys.argv[1] # .hda5 file is the first argument
  summarize_h5ad(file_path)
