#!/bin/bash

# Create directory for logs
mkdir -p submit_run_STpath_log

# Load required modules
ml fhPython/3.9.6-foss-2021b
ml TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1

# Make the script executable
chmod ug+x run_STpath.py

# Set parameters
project=Wu
task=regression
data_file=../../st2image_data/Wu_2021/data/NN_output.csv
result_path=../output/BRCA/try/
image_path=../../st2image_data/Wu_2021/output/patch_jpg/
outcome_list="Endothelial CAFs PVL B.cells T.cells Myeloid Normal.Epithelial Plasmablasts Cancer.Epithelial"
patch_id=X
base_model=ResNet50
image=o
optimizer=Adam
batch_size=64
learning_rate=0.001
dropout_rate=0
dense_layer_size=0

# Activate conda environment
source /home/zsui/miniconda3/etc/profile.d/conda.sh
conda activate tf-3.9

# Set number of threads
export OMP_NUM_THREADS=1

# Submit the job
sbatch --gpus=rtx2080ti:1 -t 10-10:00:00 -o "submit_run_STpath_log/${project}_${outcome_list}_${base_model}_${image}_${optimizer}_${batch_size}_${learning_rate}_${dropout_rate}_${dense_layer_size}.log" \
       -p campus-new --job-name=class_try \
       --wrap "python run_STpath.py \
               --project ${project} \
               --task ${task} \
               --data_file ${data_file} \
               --result_path ${result_path} \
               --image_path ${image_path} \
               --outcome_list ${outcome_list} \
               --patch_id ${patch_id} \
               --base_model ${base_model} \
               --image ${image} \
               --optimizer ${optimizer} \
               --batch_size ${batch_size} \
               --learning_rate ${learning_rate} \
               --dropout_rate ${dropout_rate} \
               --dense_layer_size ${dense_layer_size}"
