# cellxgene_handlers
Scripts for processing high dimensional data from the CZI CELLxGENE dataset 

* [The main repo](https://github.com/chanzuckerberg/single-cell-curation/tree/main) 
* [Dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md)

### Quickstart:

1) Use this command to get all available collections
```
curl -X 'GET' 'https://api.cellxgene.cziscience.com/curation/v1/collections?visibility=PUBLIC' -H 'accept: application/json' > collections.json
```

2) Extract collection IDs in collections.json with ```write_collection_ids.py```, and pass to  ```process_collections.py``` to get directory structure organized by collection ID with metadata (in JSON format) and URL information for subsequent datasets, with URL beginning in https://datasets.cellxgene.cziscience.com/ and ending in .h5ad.

3) Use ```iter_download.py``` to download the data. Optional: Use ```generate_parallel_job.py``` and ```split_for_launch.py``` for parallel download on TACC.

4) Use ```h5ad_summarizer.py``` to check out .h5ad formatted datasets.

5) When the download has finished, use ```check_h5ad.py``` to check the h5ad files in each sub-directory with their matches in each *_dataset_url.csv file to catch any files that were missed.

### Dataset structure:

```
Collections/
│
├── collection_id_0/
│   ├── collection_id_0_dataset_urls.csv
│   ├── collection_id_0_metadata.json
│   ├── dataset_id_0.h5ad
│   ├── dataset_id_1.h5ad
│   └── ...
│
├── collection_id_1/
│   ├── collection_id_1_dataset_urls.csv
│   ├── collection_id_1_metadata.json
│   ├── dataset_id_2.h5ad
│   ├── dataset_id_3.h5ad
│   └── ...
│
├── ...
│
└── collection_id_n/
    ├── collection_id_n_dataset_urls.csv
    ├── collection_id_n_metadata.json
    ├── dataset_id_x.h5ad
    ├── dataset_id_y.h5ad
    └── ...

```
Contains expression data from both humans and mice


Consolidated tables for identifying collections and datasets of interest can be found in

```./data/aggregated_metadata_tables/```
   
6) Use ```export_CLI.py``` for retrieving raw count data in csv form. Details below:

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
