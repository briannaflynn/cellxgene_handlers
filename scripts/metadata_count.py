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

        # Convert each level's Counter object to a DataFrame
        levels_dfs = []
        for level_count in self.levels_counts:
            df = pd.DataFrame(list(level_count.items()), columns=['Key', 'Count'])
            levels_dfs.append(df)

        return levels_dfs

data_processor = MetadataProcessor()
json_file_paths = data_processor.read_single_column_file_to_list('../data/metadata_files.txt')
levels_dfs = data_processor.count_items_in_json_files(json_file_paths)

for idx, df in enumerate(levels_dfs):
    print(f"DataFrame for Level {idx}:")
    print(df, "\n")
    # Save each level's DataFrame to a CSV file
    df.to_csv(f"level_{idx}.csv")
