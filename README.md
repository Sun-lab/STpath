# STpath
## Prerequisites
Before running the pipeline, ensure that the following software and libraries are installed:
   - Python 3.9.6
   - TensorFlow 2.7.1
   - Required Python libraries: numpy, pandas, argparse, ast, matplotlib, cv2, scikit-learn
   - Shell access with Slurm workload manager

## Script Breakdown
### `run_STpath.py`
This script handles the main logic for the pipeline:
1. Setup Paths: Creates necessary directories for saving results.
2. Prepare Dataframe: Adds image paths to the dataframe and splits the data into training, validation, and testing sets.
3. Validate Image Paths: Ensures all image paths are valid.
4. Create Data Generators: Generates batches of tensor image data
5. Create Model: Creates and compiles a transfer learning model.
6. Train Model: Trains the model with optional fine-tuning.
7. Evaluate Model: Evaluates the model on validation/testing datasets.
8. Plot Results: Plots learning curves and ROC curves.
9. Run Hyperparameters Tuning: Runs hyperparameter tuning and model training.
10. Main Function: Runs the training process with specified parameters.

### `submit_run_STpath.sh`
This shell script sets up the environment, loads necessary modules, activates the conda environment, and submits the job to the Slurm workload manager with the specified parameters.
1. Arguments:
   - `project`: Name of your project.
   - `task`: Type of task, either classification or regression.
   - `data_file`: Path to the CSV file containing data.
   - `result_path`: Directory where the results will be saved.
   - `image_path`: Directory containing the image patches.
   - `outcome_list`: List of columns in the CSV file, `data_file`, representing the target variables.
   - `patch_id`: Column name in `data_file` representing the image filenames or IDs.
   - `base_model`: Name of the base model for transfer learning (e.g., ResNet50).
   - `image`: Type of image normalization (o, m, v).
   - `optimizer`: Optimizer for training (e.g., Adam).
   - `batch_size`: Batch size for training.
   - `learning_rate`: Learning rate for training.
   - `dropout_rate`: Dropout rate for the dropout layer.
   - `dense_layer_size`: Number of units in the dense layer.
   - `val_split`:
   - `test_split`:
2. Example Command:
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

## Step-by-Step Guide
### Setup Environment
1. Load Required Modules: Ensure that the necessary modules are loaded. This can be done via the shell script.
2. Activate Conda Environment: Activate the conda environment where TensorFlow is installed.
### Prepare Data
1. Generate response variables.
2. Create image patches.
3. Ensure the data is in CSV format and located in the appropriate directory. The CSV file should contain image filenames or image IDs and corresponding response variables. For classification tasks, save all labels in one column; for regression tasks, save each response variable in a separate column.
### Edit Configuration and Run the Shell Script
1. Update the parameters in the shell script `submit_run_STpath.sh` according to your project.
2. Run the Shell Script. Logs will be created and saved in the submit_run_STpath_log directory.

## Output Description
After running the pipeline, you will find the following files and folders in the result path specified:
```
result_path/
├── model_project_outcome/
│   ├── model.hdf5
│   ├── fine_tuned_model.hdf5 (if fine-tuning is enabled)
│   └── ...
├── test_project_outcome/
│   ├── pred_model.csv
│   ├── eval_model.txt
│   └── ...
├── val_project_outcome/
│   ├── pred_model.csv
│   ├── eval_model.txt
│   └── ...
├── metrics_project_outcome/
│   ├── Model_Accuracy.png
│   ├── Model_Loss.png
│   └── ...
└── roc_project_outcome/
    ├── ROC_curve_val_model.png
    └── ROC_curve_test_model.png
```
Description of Files
1. model_project_outcome/: This directory contains the saved models.
   - `model.hdf5`: The primary trained model.
   - `fine_tuned_model.hdf5`: The fine-tuned model (if fine-tuning is enabled).
2. test_project_outcome/: This directory contains the results and evaluations for the test dataset.
   - `pred_model.csv`: Predictions made by the model on the test dataset.
   - `eval_model.txt`: Evaluation metrics such as classification report and confusion matrix for the test dataset.
3. val_project_outcome/: This directory contains the results and evaluations for the validation dataset.
   - `pred_model.csv`: Predictions made by the model on the validation dataset.
   - `eval_model.txt`: Evaluation metrics such as classification report and confusion matrix for the validation dataset.
4. metrics_project_outcome/: This directory contains the plots for training and validation metrics.
   - `Model_Accuracy.png`: Plot of the model accuracy over epochs.
   - `Model_Loss.png`: Plot of the model loss over epochs.
5. roc_project_outcome/: This directory contains the ROC curve plots for the validation and test datasets.
   - `ROC_curve_val_model.png`: ROC curve for the validation dataset.
   - `ROC_curve_test_model.png`: ROC curve for the test dataset.









