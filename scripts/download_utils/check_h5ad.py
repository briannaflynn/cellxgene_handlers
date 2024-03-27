import os
import csv

class H5ADFileChecker:
    def __init__(self, start_dir, suffix='_dataset_urls.csv'):
        self.start_dir = start_dir
        self.suffix = suffix

    def extract_filename_from_url(self, url):
        # This function assumes the filename can be extracted as the last part of the URL.
        # Modify this logic if the actual filename is encoded differently in the URL.
        return url.split('/')[-1]

    def check_h5ad_files(self, csv_file):
        directory = os.path.dirname(csv_file)
        missing_h5ad_urls = []

        with open(csv_file, 'r', encoding='utf-8') as file:
            urls = file.read().splitlines()
            for url in urls:
                expected_filename = self.extract_filename_from_url(url)
                # Check if the expected file exists
                if not os.path.exists(os.path.join(directory, expected_filename)):
                    missing_h5ad_urls.append(url)

        return missing_h5ad_urls

    def find_and_check_csv_files(self):
        all_missing_h5ad_urls = {}
        for root, dirs, files in os.walk(self.start_dir):
            for filename in files:
                if filename.endswith(self.suffix):
                    csv_file_path = os.path.join(root, filename)
                    missing_h5ad_urls = self.check_h5ad_files(csv_file_path)
                    if missing_h5ad_urls:
                        all_missing_h5ad_urls[csv_file_path] = missing_h5ad_urls

        return all_missing_h5ad_urls


start_directory_path = './collections/'
checker = H5ADFileChecker(start_directory_path)
missing_files_report = checker.find_and_check_csv_files()

for csv_file, missing_urls in missing_files_report.items():
    print(f"Missing files for {csv_file}:")
    for url in missing_urls:
        print(f"- {url}")
