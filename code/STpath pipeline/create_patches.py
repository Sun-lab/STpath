#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
--------------------------------------------------------------------------------
Date last modified: May 31, 2024 by Zhining Sui
Program: create_patches.py
Purpose: Create patches with exactly one filtered spot in each patch.
--------------------------------------------------------------------------------
Data Inputs:
- input_csv: Path to the CSV file with paths to required files for each WSI.
    - `sample_id`: Unique identifier for each sample/WSI.
    - `image_path`: Path to the full-resolution image file (`*.tif`).
    - `diameter`: Diameter of the spots (can be a path to `scalefactors_json.json` file, numbers, or empty).
    - `spatial_path`: Path to the spatial pixel coordinates file.
    - `barcode_path`: Path to the barcodes file, or to the filtered matrix if there is no barcodes file.
    - `output_dir`: Directory where the patches will be saved.
    - `barcode_col`: Column number for the barcode in `spatial_path`.
    - `array_row_col`: Column number for the row position in `spatial_path` (optional).
    - `array_col_col`: Column number for the column position in `spatial_path` (optional).
    - `pxl_x_col`: Column number for the x pixel coordinate in `spatial_path`.
    - `pxl_y_col`: Column number for the y pixel coordinate in `spatial_path`.

Data Outputs:
- Patches of the filtered spots saved in `output_dir`, either in JPG or TIFF format.
--------------------------------------------------------------------------------
Functions:
- locate_pixels: Overlay crosses on the image at specific coordinates.
- rgba2rgb: Convert RGBA image to RGB.
- read_sample_ids: Read sample IDs from a text file.
- find_radius: Find the radius of the spots using spatial coordinates.
- plot_diff: Plot differences between adjacent spots.
- create_patches: Create patches for each filtered spot in the specified format.
--------------------------------------------------------------------------------
Notes:
1. Adjusting Patch Radius:
    - Modify the global parameter `FACTOR` to adjust the radius of the patches relative to the radius of the spots. 
    - For 10X Visium, the radius of a spot is 27.5µm, and the distance between the centers of two spots should be 100µm. So, `FACTOR` can be at most 1.8 to avoid overlap between patches.
2. Coordinate Handling:
    - Ensure the pixel coordinates in the `spatial_path` CSV file are correctly assigned to the x and y axes. 
    - Check the generated `spots_{sample_id}.jpg` in `output_dir` to verify if the crosses match the WSI. Adjust the coordinate assignment if necessary.
    - For 10X Visium, the 5th and 6th columns in `tissue_positions_list.csv` are the pixel coordinates of the y and x axes, respectively. The coordinate system for the Python Imaging Library has the origin at the upper left corner. The x-axis (horizontal) increases from left to right, and the y-axis (vertical) increases from top to bottom. 
    - For 10X Visium, the 3rd and 4th columns in `tissue_positions_list.csv` are the row and column positions of the spot, respectively. Ensure the detected image has the correct orientation (hourglass at the upper left corner, triangle at the lower left corner, filled hexagon at the upper right corner, and open hexagon at the lower right corner). 
    - To locate spots in an image using pixel coordinates, use the `locate_pixels()` function.
    - For detailed information on fiducial alignment for spatial transcriptomics, visit:
      https://support.10xgenomics.com/spatial-gene-expression/software/visualization/latest/alignment#id_fiducials
-------------------------------------------------------------------------------
"""

import numpy as np
import tifffile as tiff
from PIL import Image
import json
import os
import pandas as pd
import gzip
import matplotlib.pyplot as plt
import seaborn as sns

# Global parameter
FACTOR = 1.1  # radius of patches/radius of spots

# function to put crosses on the image with specific coordinates
def locate_pixels(sample_id, image_path, pts, save_dir):
    """
    Overlay crosses on the image at specific coordinates.
    Args:
        sample_id (str): Unique identifier for each sample/WSI.
        image_path (str): Path to the full-resolution image file.
        pts (ndarray): Coordinates of the spots.
        save_dir (str): Directory to save the output image with crosses.
    Outputs:
        A saved image with crosses at specified coordinates.
    """

    Image.MAX_IMAGE_PIXELS = None 
    image = Image.open(image_path)
    image_np = np.array(image)
    
    width, height, channel = image_np.shape

    plt.imshow(image_np)
    # plt.plot(0, 0, "og", markersize=5)  # og:shorthand for green circle
    # plt.plot(width, height, "og", markersize=5)  # og:shorthand for green circle
    plt.scatter(pts[:, 0], pts[:, 1], marker="x", color="red", s=10)
    plt.title(sample_id)
    
    plt.savefig(f'{save_dir}/spots_{sample_id}.jpg', dpi=300)
    plt.show()
    

def rgba2rgb(rgba, background=(255, 255, 255)):
    """
    Convert RGBA image to RGB.
    Args:
        rgba (ndarray): The RGBA image.
        background (tuple): The background color to use.
    Returns:
        ndarray: The RGB image.
    """
    row, col, ch = rgba.shape
    if ch == 3:
        return rgba
    assert ch == 4, 'RGBA image has 4 channels.'
    rgb = np.zeros((row, col, 3), dtype='float32')
    a = rgba[:, :, 3] / 255.0
    for i in range(3):
        rgb[:, :, i] = rgba[:, :, i] * a + (1.0 - a) * background[i]
    return rgb.astype('uint8')

def read_sample_ids(file_path):
    """
    Read sample IDs from a text file.
    Args:
        file_path (str): Path to the text file containing sample IDs.
    Returns:
        list: List of sample IDs.
    """
    with open(file_path, 'r') as file:
        sample_ids = file.read().splitlines()
    return sample_ids

def find_radius(coords_fn = None, coords = None, diff_plot=False):
    """
    Find the radius of the spots using spatial coordinates.
    Args:
        coords_fn (str): Path to the coordinates file.
        coords (DataFrame): DataFrame containing spatial coordinates.
        diff_plot (bool): Whether to plot differences between adjacent spots.
    Returns:
        DataFrame: Updated coordinates with calculated radius.
        (Optional) Plots of differences between adjacent spots.
    """
    # Load coordinates file
    if coords is None:
        if coords_fn.endswith(".gz"):
            with gzip.open(coords_fn, 'rt') as f:
                coords = pd.read_csv(f)
        else:
            coords = pd.read_csv(coords_fn)
    
    # Rename columns and split spot column into row and col
    coords.rename(columns={coords.columns[0]: "barcode"}, inplace=True)
    coords[['array_row', 'array_col']] = coords['barcode'].str.split('x', expand=True).astype(int)
 
    adjacent_coords = pd.DataFrame(columns=["same", "spot_1", "spot_2", "dx", "dy"])  
    for axis in ['array_row', 'array_col']:
        unique_vals = coords[axis].unique()
        for val in unique_vals:
            group = coords[coords[axis] == val]
            sorted_group = group.sort_values(by=['array_col' if axis == 'array_row' else 'array_row'])
            if len(group) > 1:
                for i in range(len(sorted_group.index)-1):
                    if (sorted_group.iat[i+1, 4 if axis == 'array_row' else 3] - sorted_group.iat[i, 4 if axis == 'array_row' else 3]) == 1:
                        adjacent_coords.loc[len(adjacent_coords.index)] = [axis, 
                                                                           sorted_group.iat[i, 0], 
                                                                           sorted_group.iat[i+1, 0], 
                                                                           (sorted_group.iat[i+1, 1] - sorted_group.iat[i, 1]), 
                                                                           (sorted_group.iat[i+1, 2] - sorted_group.iat[i, 2])]
    mean_diff = adjacent_coords.groupby('same', as_index=False)[['dx', 'dy']].mean().assign(sample=coords_fn)
    var_diff = adjacent_coords.groupby('same', as_index=False)[['dx', 'dy']].var().assign(sample=coords_fn)
    
    if diff_plot:
        plot_diff(adjacent_coords, coords_fn)
    
    return coords, mean_diff, var_diff

def plot_diff(adjacent_coords, coords_fn):
    """
    Plot differences between adjacent spots.
    Args:
        adjacent_coords (DataFrame): DataFrame containing adjacent spot differences.
        coords_fn (str): Path to the coordinates file.
    Outputs:
        A plot showing differences between adjacent spots.
    """
    adjacent_coords = adjacent_coords.reset_index()
    adjacent_coords_long = pd.melt(adjacent_coords, id_vars=['index', 'same'], value_vars=['pxl_x', 'pxl_y'])

    plt.figure()
    sns.relplot(
        data=adjacent_coords_long, x="index", y="value",
        col="variable", hue="same", style="same",
        kind="scatter"
    ).set(title=coords_fn)
    plt.show()
    
    
def create_patches(samples_df, output_format = 'jpg'):
    """
    Create patches for each filtered spot in the specified format.
    Args:
        samples_df (DataFrame): DataFrame containing sample information.
        output_format (str): Format to save the patches ('jpg' or 'tiff').
    Outputs:
        Patches saved in the specified output directory.
    """
    for index, row in samples_df.iterrows():
        sample_id = row['sample_id']
        image_path = row['image_path']
        diameter = row['diameter']
        spatial_path = row['spatial_path']
        barcode_path = row['barcode_path']
        output_dir = row['output_dir']
        barcode_col = row['barcode_col']
        array_row_col = row['array_row_col']
        array_col_col = row['array_col_col']
        pxl_x_col = row['pxl_x_col']
        pxl_y_col = row['pxl_y_col']
        
        try:
            # load spatial information
            sp = pd.read_csv(spatial_path, header=None, index_col=None)
            sp.columns = [f'col_{i+1}' for i in range(sp.shape[1])]
            rename_dict = {
               f'col_{barcode_col}': 'barcode',
               f'col_{pxl_x_col}': 'pxl_x',
               f'col_{pxl_y_col}': 'pxl_y'
            }
            if pd.notna(array_row_col):
                rename_dict[f'col_{array_row_col}'] = 'array_row'
            if pd.notna(array_col_col):
                rename_dict[f'col_{array_col_col}'] = 'array_col'
            sp = sp.rename(columns=rename_dict)

            # Subset the DataFrame to only include the renamed columns
            sp = sp[list(rename_dict.values())].dropna()   

            # Ensure correct data types
            sp['barcode'] = sp['barcode'].astype(str)
            sp['pxl_x'] = sp['pxl_x'].astype(float)
            sp['pxl_y'] = sp['pxl_y'].astype(float)
            if 'array_row' in sp.columns:
                sp['array_row'] = sp['array_row'].astype(int)
            if 'array_col' in sp.columns:
                sp['array_col'] = sp['array_col'].astype(int)
            

            # Read diameter from JSON file or calculate if not available
            if pd.isna(diameter):
                sp, mean_diff, var_diff = find_radius(coords = sp) 
                # Calculate diameter: ST data has 200 µm center-to-center distance between spots, each with a diameter of 100 µm. 
                diameter = max( mean_diff[["dx", "dy"]].max(axis=1) ) / 2 
                samples_df.at[index, 'diameter'] = diameter
            elif isinstance(diameter, str) and diameter.endswith('.json'):
                try:
                    with open(diameter, 'r') as f:
                        diameter_data = json.load(f)
                        diameter = diameter_data['spot_diameter_fullres']
                except FileNotFoundError:
                    print(f'{diameter} not found.')
            else:
                diameter = diameter
            
            
            # load barcodes for the filtered matrix (in-tissue spots)
            barcodes_df = pd.read_csv(gzip.open(barcode_path, 'rt'), sep='\t', header=None, index_col=None)
            if barcodes_df.shape[1] == 1:
                barcodes = barcodes_df[0].tolist()
            else: # if there is no separate barcodes file, but only has the filtered matrix file.
                barcodes = [x for x in barcodes_df[0].tolist() if str(x) != 'nan']
            # extract spatial information of the filtered spots by matching barcodes
            sp_filtered = sp[sp["barcode"].isin(barcodes)]
            
            os.makedirs(output_dir, exist_ok=True)
            
            # check if the coordinates of spots correspond correctly with the image (image in correct orientation)
            PTS = np.array(sp_filtered[["pxl_x", "pxl_y"]])
            locate_pixels(sample_id = sample_id, image_path=image_path, pts = PTS, save_dir = output_dir)
            
            # get WSI
            Image.MAX_IMAGE_PIXELS = None
            slide = Image.open(image_path)
            slide_rgb = rgba2rgb(np.array(slide))
    
            for _, spot_row in sp_filtered.iterrows():
                # coordinates of the center of spots
                centerX, centerY = spot_row['pxl_x'], spot_row['pxl_y']
                # set the radius of the patch 
                patch_radius = diameter * FACTOR / 2
                # Crop the spot:
                x, y = centerX - patch_radius, centerY - patch_radius
                cropped_img = slide_rgb[round(y):round(y + 2 * patch_radius), round(x):round(x + 2 * patch_radius)]
                
                if output_format == "tiff":
                    tiff.imwrite(f'{output_dir}/{sample_id}_{spot_row["barcode"]}.tif', cropped_img)
                elif output_format == "jpg":
                    Image.fromarray(cropped_img).save(f'{output_dir}/{sample_id}_{spot_row["barcode"]}.jpg', "JPEG", quality=100)
                print(f'Created patch for spot: {sample_id}_{spot_row["barcode"]} in {output_format} format')
        except Exception as e:
            print(f'Error processing sample {sample_id}: {e}')
    
def main():
    INPUT_CSV = '../../data/Wu_2021/create_patches_input.csv'  
    # INPUT_CSV = '../../data/He_2020/create_patches_input.csv' 
    # NPUT_CSV = '../../data/10X/create_patches_input.csv'  
    
    samples_df = pd.read_csv(INPUT_CSV)
    
    required_columns = ['sample_id', 'image_path', 'diameter', 'spatial_path', 'barcode_path', 'output_dir',
                        'barcode_col', 'array_row_col', 'array_col_col', 'pxl_y_col', 'pxl_x_col']
    for col in required_columns:
        if col not in samples_df.columns:
            samples_df[col] = np.nan

    create_patches(samples_df, 'jpg')

    
if __name__ == '__main__':
    main()
