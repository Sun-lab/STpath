#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
--------------------------------------------------------------------------------
Date last modified: May 31, 2024 by Zhining Sui
Program: image_preprocessing.py
Purpose: Process image patches to determine if they contain a sufficient amount of tissue for further analysis.

This script uses various image processing techniques such as contrast enhancement, thresholding, and stain normalization.
The script supports multiple stain normalization methods and saves the processed images in the specified format.
--------------------------------------------------------------------------------
Data Inputs:
- dir_img_target: A target image path for stain normalization.
- dir_to_transform: A directory of images to be processed and normalized.
- output_format: The format to save normalized images (e.g., 'tiff', 'jpeg', 'jpg', 'png').
- output_path: The directory to save normalized images.
- method: The stain normalization method ('macenko', 'vahadane', 'reinhard').
- tissue_check_threshold: (Optional) The minimum tissue area ratio to consider the patch as valid.
- save_tissue_check_dir: (Optional) Directory to save the tissue check plots.

Data Outputs:
- Normalized images saved in the specified format and directory.
- (Optional) Plots showing tissue detection if tissue_check_threshold and save_tissue_check_dir are provided.
--------------------------------------------------------------------------------
Functions:
- save_image: Saves the processed image in the specified format.
- stain_normalize: Normalizes the stain of images in a directory using a specified method and optionally checks for sufficient tissue.
- is_valid_patch: Checks if an image patch contains sufficient tissue area for further analysis.
"""

import os
import re
import cv2
import numpy as np
import torch
from torchvision.transforms import ToTensor
from torchvision.transforms.functional import convert_image_dtype
from torch_staintools.normalizer import NormalizerBuilder
from torch_staintools.augmentor import AugmentorBuilder
from PIL import Image
import tifffile as tiff
import matplotlib.pyplot as plt
from skimage import filters, exposure, morphology


def save_image(image, output_format, output_path, image_id):
    """
    Save the processed image in the specified format.

    Args:
        image (ndarray): The image to save.
        output_format (str): The format to save the image ('tiff', 'jpeg', 'jpg', 'png').
        output_path (str): The path to save the image.
        image_id (str): The identifier for the image.
    """
    if output_format == "tiff":
        tiff.imwrite(os.path.join(output_path, f'{image_id}.tif'), image)
    elif output_format in ["jpeg", "jpg"]:
        Image.fromarray(image).save(os.path.join(output_path, f'{image_id}.jpg'), "JPEG", quality=100)
    elif output_format == "png":
        Image.fromarray(image).save(os.path.join(output_path, f'{image_id}.png'), "PNG", quality=100)
    else:
        raise ValueError(f"Unsupported output format: {output_format}")

def stain_normalize(dir_img_target, dir_to_transform, output_format, output_path, method, tissue_check_threshold=0, save_tissue_check_dir=None):
    """
    Perform stain normalization on images in a directory.

    Args:
        dir_img_target (str): The target image path for normalization.
        dir_to_transform (str): The directory of images to normalize.
        output_format (str): The format to save normalized images.
        output_path (str): The path to save normalized images.
        method (str): The normalization method ('macenko', 'vahadane', 'reinhard').
        tissue_check_threshold (float): The minimum tissue area ratio to consider the patch as valid.
        save_tissue_check_dir (str): Directory to save the tissue check plots.
    """
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    img_format = '.' + dir_img_target.split("/")[-1].split(".")[-1]
    
    seed = 0
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    # cpu or gpu
    device = torch.device("cpu")
    
    # shape: Height (H) x Width (W) x Channel (C, for RGB C=3)
    target = cv2.cvtColor(cv2.imread(dir_img_target), cv2.COLOR_BGR2RGB)
    # shape: Batch x Channel x Height x Width (BCHW)
    target_tensor = ToTensor()(target).unsqueeze(0).to(device)
    
    # ######## Normalization
    # Create the normalizer
    if method.lower() == 'macenko':
        normalizer = NormalizerBuilder.build('macenko', concentration_method='ista')
    elif method.lower() == 'vahadane':
        normalizer = NormalizerBuilder.build('vahadane', concentration_method='ista')
    elif method.lower() == 'reinhard':
        normalizer = NormalizerBuilder.build('reinhard', concentration_method='ista')
    else:
        raise ValueError(f"Unsupported normalization method: {method}")

    # Move the normalizer to the device (CPU or GPU)
    normalizer = normalizer.to(device)
    # Fit. For macenko and vahadane this step will compute the stain matrix and concentration
    normalizer.fit(target_tensor)

    if os.path.isdir(dir_to_transform):
        for img in os.listdir(dir_to_transform):
            if img.endswith(img_format) and not os.path.isfile(output_path + img):
                to_transform_id = re.sub(img_format, '', img)
                # shape: Height (H) x Width (W) x Channel (C, for RGB C=3)
                to_transform = cv2.cvtColor(cv2.imread(os.path.join(dir_to_transform, img)), cv2.COLOR_BGR2RGB)
                if is_valid_patch(to_transform, save_dir=save_tissue_check_dir, patch_id=to_transform_id, tissue_threshold=tissue_check_threshold):
                    # transform
                    to_transform_tensor = ToTensor()(to_transform).unsqueeze(0).to(device)
                    # BCHW - scaled to [0, 1] torch.float32
                    norm = normalizer(to_transform_tensor)
                    norm_np = convert_image_dtype(norm, torch.uint8).squeeze().detach().cpu().permute(1, 2, 0).numpy()                    
                    save_image(norm_np, output_format, output_path, to_transform_id)
                    print(f"The image {img} has been normalized and saved in the destination directory.")
            else:
                print(f"The image {img} has already been normalized and saved in the destination directory.")
    elif os.path.isfile(dir_to_transform):
        to_transform_id = os.path.basename(dir_to_transform).split(".")[0]
        if not os.path.isfile(output_path + to_transform_id + ".jpg"):
            to_transform = cv2.cvtColor(cv2.imread(dir_to_transform), cv2.COLOR_BGR2RGB)
            if is_valid_patch(to_transform, save_dir=save_tissue_check_dir, patch_id=to_transform_id, tissue_threshold=tissue_check_threshold):
                # transform
                to_transform_tensor = ToTensor()(to_transform).unsqueeze(0).to(device)
                # BCHW - scaled to [0, 1] torch.float32
                norm = normalizer(to_transform_tensor)
                norm_np = convert_image_dtype(norm, torch.uint8).squeeze().detach().cpu().permute(1, 2, 0).numpy()                    
                save_image(norm_np, output_format, output_path, to_transform_id)
        print(f"The image {os.path.basename(dir_to_transform)} has been normalized and saved in the destination directory.")


def is_valid_patch(patch, tissue_threshold=0.5, contrast_method='histogram', threshold_method='hysteresis', hist_params=None, hyst_params=None, clean_params=None, display=False, save_dir=None, patch_id="patch"):
    """
    Check if the patch is valid by analyzing the tissue area.

    Args:
        patch (ndarray): The input image patch.
        tissue_threshold (float): The minimum tissue area ratio to consider the patch as valid.
        contrast_method (str): The contrast enhancement method.
        threshold_method (str): The thresholding method.
        hist_params (dict): Parameters for the histogram equalization method.
        hyst_params (dict): Parameters for the hysteresis thresholding method.
        clean_params (dict): Parameters for the function remove_small_objects().
        display (bool): Whether to display the intermediate results.
        save_dir (str): Directory to save the plots.
        patch_id (str): ID of the patch to be used in the saved figure filename.

    Returns:
        bool: True if the tissue area ratio is above the threshold, False otherwise.
    """
    
    # Default parameters for contrast methods
    if hist_params is None:
        hist_params = {'nbins': 256, 'clip_limit': 0.05}
    if hyst_params is None:
        hyst_params = {'low': 140, 'high': 230}
    if clean_params is None:
        clean_params = {'min_size': 300}

    # Convert the patch to grayscale
    gray_patch = cv2.cvtColor(patch, cv2.COLOR_RGB2GRAY)

    # Increase contrast in the grayscale patch if specified
    if contrast_method == 'histogram':
        contrast_patch = exposure.equalize_hist(gray_patch, nbins=hist_params['nbins'])
        contrast_patch = (contrast_patch * 255).astype(np.uint8)  # Convert to 8-bit
    elif contrast_method == 'adaptive':
        contrast_patch = exposure.equalize_adapthist(gray_patch, nbins=hist_params['nbins'], clip_limit=hist_params['clip_limit'])
        contrast_patch = (contrast_patch * 255).astype(np.uint8)  # Convert to 8-bit
    elif contrast_method == 'none':
        contrast_patch = gray_patch
    else:
        raise ValueError("Invalid contrast method. Choose 'none', 'histogram', or 'adaptive'.")

    # Obtain the complement of the increased contrast grayscale image
    complement_patch = cv2.bitwise_not(contrast_patch)
    
    # Apply the chosen thresholding method
    if threshold_method == 'hysteresis':
        hyst = filters.apply_hysteresis_threshold(complement_patch, low=hyst_params['low'], high=hyst_params['high'])
        cleaned_threshold = morphology.remove_small_objects(hyst, min_size=clean_params['min_size'])
        # cleaned_threshold = morphology.binary_fill_holes(cleaned_threshold)
    elif threshold_method == 'otsu':
        otsu_thresh_value = filters.threshold_otsu(complement_patch)
        otsu = (complement_patch > otsu_thresh_value)
        cleaned_threshold = morphology.remove_small_objects(otsu, min_size=clean_params['min_size'])
        # cleaned_threshold = binary_fill_holes(cleaned_threshold)
    else:
        raise ValueError("Invalid threshold method. Choose 'hysteresis' or 'otsu'.")

    # Calculate tissue area: count pixels in filled regions that are not black in the complement image
    tissue_area = np.sum(cleaned_threshold)
    total_area = complement_patch.size
    tissue_ratio = tissue_area / total_area

    # Display and save the original image, grayscale image, complement image, edges, filled regions, and tissue area if requested
    if display or save_dir:
        fig, axs = plt.subplots(1, 6, figsize=(25, 5))
        
        # Display original image
        axs[0].imshow(patch)
        axs[0].set_title('Original Image')
        axs[0].axis('off')
        
        # Display grayscale image
        axs[1].imshow(gray_patch, cmap='gray')
        axs[1].set_title('Grayscale Image')
        axs[1].axis('off')
        
        # Display contrast enhanced image if applicable
        axs[2].imshow(contrast_patch, cmap='gray')
        if contrast_method == 'none':
            axs[2].set_title('Original Grayscale Image')
        else:
            axs[2].set_title('Contrast Enhanced Image')
        axs[2].axis('off')
        
        # Display complement image
        axs[3].imshow(complement_patch, cmap='gray')
        axs[3].set_title('Complement Image')
        axs[3].axis('off')
        
        # Display thresholded image
        axs[4].imshow(cleaned_threshold, cmap='gray')
        if threshold_method == 'hysteresis':
            axs[4].set_title('Hysteresis Thresholded')
        else:
            axs[4].set_title('Otsu Thresholded')
        axs[4].axis('off')
        
        # Display tissue area
        tissue_overlay = np.zeros_like(patch)
        tissue_overlay[cleaned_threshold] = [255, 0, 0]  # Highlight tissue area in red
        combined_image = cv2.addWeighted(patch, 0.7, tissue_overlay, 0.3, 0)
        axs[5].imshow(combined_image)
        axs[5].set_title('Tissue Area')
        axs[5].axis('off')
        
        if display:
            plt.show()
        
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            fig.savefig(os.path.join(save_dir, f'{patch_id}_tissue_detected.png'))
            plt.close(fig)
            
    # Determine if the patch is valid
    return tissue_ratio >= tissue_threshold




# Example usage
ORIGINAL_IMG_DIR = '/Users/zhiningsui/GitHub/STpath/output/10X/patch/'
NORMALIZED_IMG_DIR = '/Users/zhiningsui/GitHub/STpath/output/10X/patch_normalized/'
TISSUE_IMG_DIR = '/Users/zhiningsui/GitHub/STpath/output/10X/patch_tissue_detected/'

for method in ['macenko', 'vahadane', 'reinhard']:
    stain_normalize(dir_img_target=ORIGINAL_IMG_DIR + '10X_FFPE_GGAACCGTGTAAATTG-1.jpg',
                    dir_to_transform=ORIGINAL_IMG_DIR,
                    output_format='jpg', output_path=NORMALIZED_IMG_DIR + method +' /', 
                    method=method, tissue_check_threshold = 0.2, save_tissue_check_dir = TISSUE_IMG_DIR)





