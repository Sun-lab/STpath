#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 16:17:55 2025

@author: zhiningsui
"""

from huggingface_hub import snapshot_download
from huggingface_hub import hf_hub_download
import pandas as pd
import re
import os

# reading csv file 
meta = pd.read_csv("/Users/zhiningsui/GitHub/STpath/data/STImage-1K4M/create_patches_input_brca.csv")

# Repository information
repo_id = "jiawennnn/STimage-1K4M"  # Replace with the correct repo ID
repo_type = "dataset"  # Specify that it's a dataset

# List of file paths within the repository
file_paths = meta["count_matrix_dir"].tolist()
file_paths = [re.sub(r"^.*STimage-1K4M/", "", x) for x in file_paths]

# Local directory to save downloaded files
local_dir = "downloaded_files"
os.makedirs(local_dir, exist_ok=True)

# Download files
for file_path in file_paths:
    try:
        local_file = hf_hub_download(
            repo_id=repo_id,
            filename=file_path,
            repo_type=repo_type,
            local_dir=local_dir
        )
        print(f"Downloaded: {file_path} to {local_file}")
    except Exception as e:
        print(f"Failed to download {file_path}: {e}")


    
    
    