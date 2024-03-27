
## Metadata structure notes

At the top level of each metadata JSON file, there's a key named datasets which contains a list of dictionaries. Each dictionary within this list includes several keys:

dataset_id: A string value that uniquely identifies the dataset.

developmental_stage: A list of dictionaries, each containing a key named label.

disease: A list of dictionaries, each containing a key named label.

organism: A list of dictionaries, each containing a key named label.

raw_data_location: A string value indicating the location of the raw data.

See metadata_organization.txt for some additional free hand notes on file structure