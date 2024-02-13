#!/usr/bin/python
import pandas as pd
from request_collection import fetch_and_save_dataset_urls

input_file_path = 'collection_ids.txt'
bad_collecs = []
with open(input_file_path, 'r') as file:
    for line in file:
        collection_id = line.strip()  
        if collection_id:  
            print(f"Processing collection_id: {collection_id}")
            try:
                fetch_and_save_dataset_urls(collection_id, create_directory=True)
            except Exception as e:
                print(f"An error occurred while processing {collection_id}: {e}")
                bad_collecs.append(collection_id)

# log any collection IDs that had issues while processing
df = pd.DataFrame()
df['bad_collection_ids'] = bad_collecs
df.to_csv('bad_collection_ids.csv', index=False)
