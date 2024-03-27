import os
import requests

class DatasetDownloader:
    def __init__(self, start_dir, suffix='_dataset_urls.csv'):
        self.start_dir = start_dir
        self.suffix = suffix

    def download_file(self, url, directory):
        try:
            response = requests.get(url, allow_redirects=True)
            # Try to get filename from content-disposition header
            filename = url.split('/')[-1]  # Default filename
            if 'content-disposition' in response.headers:
                content_disposition = response.headers['content-disposition']
                filename_parts = content_disposition.split('filename=')
                if len(filename_parts) > 1:
                    filename = filename_parts[1]
                    if filename[0] == '"' or filename[0] == "'":
                        filename = filename[1:-1]
            # Ensure directory exists
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, filename)
            # Save the file
            with open(filepath, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded {filename} to {directory}")
        except Exception as e:
            print(f"Error downloading {url}: {e}")

    def process_csv_file(self, csv_file):
        directory = os.path.dirname(csv_file)
        with open(csv_file, mode='r', encoding='utf-8') as file:
            urls = file.read().splitlines()
            for url in urls:
                self.download_file(url, directory)

    def find_and_process_csv_files(self):
        for root, dirs, files in os.walk(self.start_dir):
            for filename in files:
                if filename.endswith(self.suffix):
                    csv_file_path = os.path.join(root, filename)
                    print(f"Processing file: {csv_file_path}")
                    self.process_csv_file(csv_file_path)