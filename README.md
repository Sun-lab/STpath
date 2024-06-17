# STpath
## Scripts
### `run_STpath.py`
This script handles the main logic for the pipeline:
1. Setup Paths: Creates necessary directories for saving results.
2. Prepare Dataframe: Adds image paths to the dataframe and splits the data into training, validation, and testing sets.
3. Validate Image Paths: Ensures all image paths are valid.
4. Create Data Generators: Generates batches of tensor image data
5. Create Model: Creates and compiles a transfer learning model.
6. Train Model: Trains the model with optional fine-tuning.
7. Predict New Data: Makes predictions on the test dataset and saves results.
8. Plot Results: Plots learning curves and ROC curves.
9. Run Hyperparameters Tuning: Runs hyperparameter tuning and model training.
10. Main Function: Runs the training process with specified parameters.
    
To run the script, use the provided shell script `submit_run_STpath.sh`. This shell script sets up the environment, loads necessary modules, activates the conda environment, and submits the job to the Slurm workload manager with the specified parameters.
- Arguments:
   - `project`: Name of your project.
   - `task`: Type of task, either "classification" or "regression".
   - `data_file`: Path to the CSV file containing data.
   - `result_path`: Directory where the results will be saved.
   - `image_path`: Directory containing the image patches.
   - `outcome_list`: List of columns in the CSV file, `data_file`, representing the target variables.
   - `patch_id`: Column name in `data_file` representing the image filenames or IDs.
   - `base_model`: Name of the base model for transfer learning (default: ResNet50).
   - `image`: Type of image normalization (o = original, m = macenko, v = vahadane).
   - `optimizer`: Optimizer for training (default: Adam).
   - `batch_size`: Batch size for training.
   - `learning_rate`: Learning rate for training.
   - `dropout_rate`: Dropout rate for the dropout layer.
   - `dense_layer_size`: Number of units in the dense layer.
   - `val_split`: Split ratio of the validation dataset (default: 0.15).
   - `test_split`: Split ratio of the testing dataset (default: 0.15).
   - `seed`: Random seed for reproducibility (default: 17).
   - `image_resize`: Size of image input to the model (default: 224).
   - `num_epoch`: Maximum number of epochs (default: 500).
   - `patience`: Early stopping patience (default: 20).
- Example command for job submission:
```
sbatch --gpus=rtx2080ti:1 -t 10-10:00:00 -o "submit_run_STpath_log/${project}_${base_model}_${outcome_list}_${image}_${optimizer}_${batch_size}_${learning_rate}_${dropout_rate}_${dense_layer_size}.log" \
       -p campus-new --job-name=class_try \
       --wrap "python run_STpath.py \
               --project He \
               --task regression \
               --data_file BRCA_proportions.csv \
               --result_path ../output/BRCA/try/ \
               --image_path ../../st2image_data/BRCA/output/patch_jpg_all/ \
               --outcome_list invasive.cancer stroma lymphocytes others \
               --patch_id X \
               --base_model ResNet50 \
               --image o \
               --optimizer Adam \
               --batch_size 64 \
               --learning_rate 0.001 \
               --dropout_rate 0 \
               --dense_layer_size 0"
```
### `create_patches.py`
This script creates patches from the full-resolution WSIs based on the ST data such that there is one spot in each patch. 
### `image_preprocess.py`
This script evaluates the quality of patches by calculating the proportion of tissue in the patch and performs stain normalization.
### `clustering_preprocess.R`
This script handles the preprocessing of the ST count matrices to prepare for clustering. It performs filtering and creates Seurat objects from the count matrices of all samples, which are ready for clustering (`clustering.R`)
### `clustering.R`
This script handles the normalization, integration, clustering, and identifying markers for a pre-processed Seurat object obtained by `clustering_preprocess.R`.
## Step-by-Step Guide
### Setup Environment
1. Load Required Modules: Ensure that the necessary modules are loaded. This can be done via the shell script.
2. Activate Conda Environment: Activate the conda environment where TensorFlow is installed.
### Prepare Data
1. Generate response variables. Save the response variable to the `data_file` CSV file. 
  - Run CARD (`.R`) to obtain the cell type proportions for regression tasks.
  - Run Seurat (`clustering_preprocess.R` and `clustering.R`) to obtain the clusters for classification tasks.
  - Use other types of response variables of interest.
2. Create image patches. Save the filename of each patch in the `data_file` CSV file if the filenames are not the same as the patch IDs. 
  - Run `create_patches.py` to generate image patches. 
4. Ensure the data is in CSV format and located in the appropriate directory. The CSV file should contain patch IDs or image patch filenames (if not the same as patch ID) and corresponding response variables. For classification tasks, save all labels in one column; for regression tasks, save each response variable in a separate column.
### Edit Configuration and Run the Shell Script
1. Update the parameters in the shell script `submit_run_STpath.sh` according to your project.
2. Run the Shell Script. Logs will be created and saved in the submit_run_STpath_log directory.

## Output Description
After running the pipeline, you will find the following files and folders in the result path specified:
```
result_path/
├── training_{project}_{outcome_list}.csv
├── testing_{project}_{outcome_list}.csv
├── validation_{project}_{outcome_list}.csv
├── model_{project}_{outcome_list}/
│   ├── {model_name}.hdf5
│   ├── fine_tuned_{model_name}.hdf5 (if fine-tuning is enabled)
│   └── ...
├── curve_{project}_{outcome_list}/
│   ├── 
│   └── ...
├── eval_{project}_{outcome_list}/
│   ├── test_pred_{model_name}.csv
│   ├── test_metrics_{model_name}.csv
|   ├── metrics_history_{model_name}.csv
|   ├── report_{model_name}.txt (if classification)
│   └── ...
└── roc_{project}_{outcome_list}/ (if classification)
    ├── ROC_curve_test_{model_name}.png
    └── ...
```
Description of Files
1. `training_{project}_{outcome_list}.csv`, `testing_{project}_{outcome_list}.csv`, `validation_{project}_{outcome_list}.csv`: training, testing, and validation datasets used in the job.
2.  model_{project}_{outcome_list}/: This directory contains the saved models.
   - `{model_name}.hdf5`: The trained model.
   - `fine_tuned_{model_name}.hdf5`: The fine-tuned model (if fine-tuning is enabled).
3. curve_{project}_{outcome_list}/: This directory contains the plots for the learning curve for training and validation.
   - `{model_name}.png`: Plots of the learning curve.
4. eval_{project}_{outcome_list}/: This directory contains the results and evaluations for the testing dataset.
   - `metrics_history_{model_name}.csv`: History of model metrics and model loss during training.
   - `test_pred_{model_name}.csv`: Predictions (scores and predicted labels for classification tasks; predicted response variables for regression tasks) made by the model on the testing dataset.
   - `test_metrics_{model_name}.csv`: Evaluation metrics (model metrics and model loss) for the testing dataset.
   - `report_{model_name}.txt`: Classification report and confusion matrix of the testing dataset (if classification). 
5. roc_{project}_{outcome_list}/: This directory contains the ROC curve plots for the test datasets.
   - `ROC_curve_test_{model_name}.png`: ROC curve for the test dataset.










