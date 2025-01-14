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


    
import os
from synapseclient import Synapse
import pandas as pd
syn = Synapse()
syn.login(authToken='eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTczNjgwNzg2NCwiaWF0IjoxNzM2ODA3ODY0LCJqdGkiOiIxNTQyNiIsInN1YiI6IjM1Mjc4OTYifQ.i4lW4HRwD0uIAX_oV0XTnwI0oSJk75FiiqbqtJdggU922FYWoGsctkHyYc_lGq9wjdVhhkwVRuV_mNp0rmg9ecGd6VfzlMDQZm6CBAlB0nowHcy8B5TERRgFY2rq8RoQ8TPV-89iCVYsxiyRqeUjFPOrOoL4lY_FkS-OkCiUDVZ-fbrI5CIhji5JOul7HPuDbfAeB0yHqEklat-a1L13N0EzOpFR5F3hf_qFFU1llubtWG2D33les3bZAcekTHMv8x-C3gKz0zJzxd5tW-VapGiUw_yqZgBBsPG36Pv5ab1PWfbks4cIp1h437baoVG4yE-o6LnXeb3kDn9Fy11DOQ')

# Load the full list of Synapse IDs from the CSV file
file_path = "/Users/zhiningsui/GitHub/STpath/code/STpath pipeline/Unsuccessful_Downloads.csv"
synapse_ids = pd.read_csv(file_path)["SynapseID"].tolist()

download_directory = "/Users/zhiningsui/GitHub/Mo_2024"  # Replace with your desired path
os.makedirs(download_directory, exist_ok=True)  # Create the directory if it doesn't exist

# Track unsuccessful downloads
unsuccessful_downloads = []

# Attempt to download each file
for syn_id in synapse_ids:
    try:
        print(f"Downloading {syn_id} to {download_directory}...")
        syn.get(syn_id, downloadLocation=download_directory)
        print(f"Successfully downloaded {syn_id}.")
    except Exception as e:
        print(f"Error downloading {syn_id}: {e}")
        unsuccessful_downloads.append({"SynapseID": syn_id, "Error": str(e)})

# Save the list of unsuccessful downloads to a CSV file
if unsuccessful_downloads:
    error_file_path = "/Users/zhiningsui/GitHub/STpath/code/STpath pipeline/Unsuccessful_Downloads.csv"
    pd.DataFrame(unsuccessful_downloads).to_csv(error_file_path, index=False)
    print(f"List of unsuccessful downloads saved to {error_file_path}.")
else:
    print("All files downloaded successfully.")




















