#!/usr/bin/bash

while IFS= read -r FILE_PATH; do python h5ad_summarizer.py "$FILE_PATH"; done < ../data/h5ad_files.txt
