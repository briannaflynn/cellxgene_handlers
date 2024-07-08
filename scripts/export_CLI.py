import argparse
from export_raw_data import *
from cxg_logger import *
import sys

docstring = """
Command line interface for raw count matrix export using export_raw_data.py

The main arguments are csv_request, data_path, and export_path

    csv_request: path to a request file with requested collection ID, dataset ID, etc, formatted like so - 

        collection_id,dataset_id,dataset_version_id,development_stage,disease,organism,raw_data_location,cell_type
        283d65eb-dd53-496d-adb7-7570c7caa443,ff7d15fa-f4b6-4a0e-992e-fd0c9d088ded,51e05270-1f00-4527-b73a-f770a5a14f62,50-year-old human stage,normal,Homo sapiens,X,neuron
        283d65eb-dd53-496d-adb7-7570c7caa443,fe1a73ab-a203-45fd-84e9-0f7fd19efcbd,4e124ecc-7885-465c-bab9-4e94d9d40b6a,50-year-old human stage,normal,Homo sapiens,X,neuron
        283d65eb-dd53-496d-adb7-7570c7caa443,fbf173f9-f809-4d84-9b65-ae205d35b523,5a52f557-aeaf-4fc9-a4a9-7a2eaca5cd3c,50-year-old human stage,normal,Homo sapiens,X,neuron
        283d65eb-dd53-496d-adb7-7570c7caa443,fa554686-fc07-44dd-b2de-b726d82d26ec,6606e9aa-e4c4-4522-9e38-8128aa415b15,50-year-old human stage,normal,Homo sapiens,X,neuron

        Absolutely necessary fields are: collection_id, dataset_version_id, and raw_data_location
    
    data_path: path to the directory containing the collection_id directories, which contain the dataset .h5ad files
        
        For example, the data path on the Marcotte pod is: /stor/scratch/External/CELLxGENE/collections/

    export_path: path to a directory where the csv files will be exported to

        This is a location of your choosing, for example - /stor/scratch/External/CELLxGENE/my_data_export/

The optional arguments are exclude IDs adn gene_filter_list

    exclude IDs: A list of dataset version IDs to exclude from the analysis

    gene_filter_list: Path to a csv filter containing the desired genes to keep, when provided all other genes that do not match are excluded

        Default is the Humap 2.0 genes listed in ../data/all_ensembl_ids_df.csv

        Format as a single column, like so -
        
        ENSEMBL_ID
        ENSG00000113552
        ENSG00000233917ß
        ENSG00000254709
        ENSG00000010017
        ENSG00000113810
        ENSG00000102103
        ENSG00000241553
        ENSG00000111229
        ENSG00000130429

"""

def main():
    parser = argparse.ArgumentParser(description=docstring, formatter_class = argparse.RawTextHelpFormatter, argument_default=argparse.SUPPRESS)
    parser.add_argument('csv_request', type=str, help="Path to the CSV request file")
    parser.add_argument('data_path', type=str, help="Top level directory containing collections/datasets")
    parser.add_argument('export_path', type=str, help="Directory where the files will be exported")
    parser.add_argument('--exclude_IDs', type=str, nargs='*', default=None, help="List of dataset version IDs to exclude")
    parser.add_argument('--gene_filter_list', type=str, default=None, help="Path to gene filter list CSV file")

    # prints description docstring and help when user provides no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    df = pd.read_csv(args.csv_request)

    # supply a csv file list of genes
    if args.gene_filter_list is not None:
        gene_filter_path = args.gene_filter_list if args.gene_filter_list else '../data/all_ensembl_ids_df.csv'
        gene_filter_df = pd.read_csv(gene_filter_path)
        keep_list = gene_filter_df['ENSEMBL_ID'].to_list()
    else:
        keep_list = None

    print('#' * 25, 'STARTING EXPORT', '#' * 25)

    date = str(dt.datetime.now())
    date = date.replace(' ', '_')

    logger = setup_logger('export', f'raw_data_export_{date}.log')

    with StreamToLogger(logger, logging.INFO):
        try:
            export_raw_counts_to_csv(df, base_path_string=args.data_path, out_path=args.export_path, exclude_IDs=args.exclude_IDs, gene_filter_list=keep_list)
        except BaseException as error:
            print(f'BASE EXCEPTION: {error}, ID was NOT processed, most likely a memory issue')
            pass

if __name__ == "__main__":
    main()