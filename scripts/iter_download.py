import pandas as pd
from download_utility import DatasetDownloader
import os

input_file_path = '../data/collection_ids.txt'
bad_collecs = []
with open(input_file_path, 'r') as file:
    for line in file:
        collection_id = line.strip()  
        if collection_id:  
            print(f"Processing collection_id: {collection_id}")
            directory_path = os.path.join(os.getcwd(), 'collections', collection_id)
            print(directory_path)
            try:
                downloader = DatasetDownloader(directory_path)
                downloader.find_and_process_csv_files()
            except Exception as e:
                print(f"An error occurred while processing {collection_id}: {e}")
                bad_collecs.append(collection_id)

print(bad_collecs)