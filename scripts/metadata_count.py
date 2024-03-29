import json
import pandas as pd
from collections import Counter

# class for managing hierarchical evaluation of metadata JSON files
class MetadataProcessor:
    def __init__(self):
        self.levels_counts = []

    def read_single_column_file_to_list(self, filename):
        with open(filename, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
        return lines

    def process_dict(self, data, level=0):
        if level >= len(self.levels_counts):
            self.levels_counts.append(Counter())
        for key, value in data.items():
            self.levels_counts[level].update([key])
            if isinstance(value, dict):
                self.process_dict(value, level + 1)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, (str, int, float)):
                        self.levels_counts[level].update([item])
                    elif isinstance(item, dict):
                        self.process_dict(item, level + 1)
            else:
                self.levels_counts[level].update([value])

    def count_items_in_json_files(self, json_file_paths):
        self.levels_counts = []  # Reset for each call to ensure fresh counts
        for file_path in json_file_paths:
            with open(file_path, 'r') as file:
                data = json.load(file)
                self.process_dict(data)

        levels_dfs = []
        for level_count in self.levels_counts:
            df = pd.DataFrame(list(level_count.items()), columns=['Key', 'Count'])
            levels_dfs.append(df)

        return levels_dfs

    def aggregate_datasets_info(self, json_file_paths, cell_type = True, verbose = True, output_csv_path='../data/aggregated_metadata_json_with_celltype.csv'):

        def fill_missing_values(df, column_name:str):

            last_valid_value = None
            for i in range(len(df)):
                if pd.notnull(df.at[i, column_name]) and df.at[i, column_name].strip():
                    last_valid_value = df.at[i, column_name]
                elif last_valid_value is not None:
                    df.at[i, column_name] = last_valid_value

                if i < len(df) - 1 and df.at[i, 'dataset_id'] != df.at[i + 1, 'dataset_id']:
                    last_valid_value = None

            return df
       
        aggregated_data = []
        columns_to_fill = ['development_stage', 'disease', 'organism']
        if cell_type:
            columns_to_fill.append('cell_type')

        for file_path in json_file_paths:
            with open(file_path, 'r') as file:
                data = json.load(file)
                collection_id = data['collection_id']
                print(f"Processing {collection_id}")
                print(file_path)
                for dataset in data['datasets']:
                    
                    dataset_id = dataset['dataset_id']

                    raw_data_location = dataset['raw_data_location']

                    if cell_type:
                        developmental_stages = dataset['development_stage']
                        diseases = dataset['disease']
                        organisms = dataset['organism']
                        cells = dataset['cell_type']

                        max_length = max(len(developmental_stages), len(diseases), len(organisms), len(cells))
                        
                        for i in range(max_length):
                            dev_stage_label = developmental_stages[i].get('label', '') if i < len(developmental_stages) else ''
                            disease_label = diseases[i].get('label', '') if i < len(diseases) else ''
                            organism_label = organisms[i].get('label', '') if i < len(organisms) else ''
                            cell_label = cells[i].get('label', '') if i < len(cells) else ''

                            aggregated_data.append({
                                'collection_id':collection_id,
                                'dataset_id':dataset_id,
                                'development_stage':dev_stage_label,
                                'disease':disease_label,
                                'organism':organism_label,
                                'raw_data_location':raw_data_location,
                                'cell_type':cell_label
                            })
                    else:
                        developmental_stages = dataset['development_stage']
                        diseases = dataset['disease']
                        organisms = dataset['organism']

                        max_length = max(len(developmental_stages), len(diseases), len(organisms))
                        
                        for i in range(max_length):
                            dev_stage_label = developmental_stages[i].get('label', '') if i < len(developmental_stages) else ''
                            disease_label = diseases[i].get('label', '') if i < len(diseases) else ''
                            organism_label = organisms[i].get('label', '') if i < len(organisms) else ''

                            aggregated_data.append({
                                'collection_id':collection_id,
                                'dataset_id':dataset_id,
                                'development_stage':dev_stage_label,
                                'disease':disease_label,
                                'organism':organism_label,
                                'raw_data_location':raw_data_location
                            })

        df = pd.DataFrame(aggregated_data)

        for column in columns_to_fill:
            df = fill_missing_values(df, column)

        if verbose:
            print(df)

        df.to_csv(output_csv_path, index=False)

        return df

data_processor = MetadataProcessor()
json_file_paths = data_processor.read_single_column_file_to_list('../data/metadata_files.txt')

data_no_celltypes = data_processor.aggregate_datasets_info(json_file_paths, cell_type=False, output_csv_path='aggregated_metadata_json.csv')
df = data_processor.aggregate_datasets_info(json_file_paths)
