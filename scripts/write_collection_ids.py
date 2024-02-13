#!/usr/bin/python

import json

json_file_path = 'collections.json'
output_file_path = 'collection_ids.txt'

# Open and read the JSON file
with open(json_file_path, 'r') as json_file:
    data = json.load(json_file)

# Open the output file in write mode
with open(output_file_path, 'w') as output_file:
    # Iterate through each dictionary in the list
    for item in data:
        # Extract the 'collection_id' value and write it to the output file
        collection_id = item['collection_id']
        output_file.write(collection_id + '\n')

print(f"'collection_id' values have been written to {output_file_path}")
