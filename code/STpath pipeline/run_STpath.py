#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
import os
import itertools
import matplotlib.pyplot as plt
import cv2


from tensorflow.keras import layers, models, optimizers, utils, applications, preprocessing, callbacks
from sklearn import metrics as skmetrics
from sklearn.utils import compute_class_weight


# Constants
POOLING_AVERAGE = 'avg'
BATCH_SIZE_VALIDATION = 1
NUM_FINE_TUNE_EPOCHS = 3

# Paths (global, set in main)
EVAL_PATH = PLOT_PATH = MODEL_PATH = ROC_PATH = ""


def setup_paths(result_path, project, outcome_list, task):
    """
    Set up the directory paths for saving results.
    
    Args:
        result_path (str): Base path for saving results.
        project (str): Project name.
        outcome (str): Outcome name.
        task (str): Task type (classification or regression)
    
    Returns:
        None
    """
    outcome_str = '_'.join(outcome_list) if isinstance(outcome_list, list) else outcome_list
    
    global EVAL_PATH, PLOT_PATH, MODEL_PATH, ROC_PATH
    EVAL_PATH = os.path.join(result_path, f"eval_{project}_{outcome_str}/")
    PLOT_PATH = os.path.join(result_path, f"curve_{project}_{outcome_str}/")
    MODEL_PATH = os.path.join(result_path, f"model_{project}_{outcome_str}/")
    os.makedirs(EVAL_PATH, exist_ok=True)
    os.makedirs(PLOT_PATH, exist_ok=True)
    os.makedirs(MODEL_PATH, exist_ok=True)
    if task == "classification":
        ROC_PATH = os.path.join(result_path, f"roc_{project}_{outcome_str}/")
        os.makedirs(ROC_PATH, exist_ok=True)
    print(f'Results are saved in: {MODEL_PATH}, {EVAL_PATH}, {PLOT_PATH}')
    

# Add image paths to dataframe and split into training, validation, and testing sets
def prepare_dataframe(df, filename, path, train_split, val_split, test_split):
    """
    Adds a column for image paths and splits the dataframe into training, validation, and testing sets.
    
    Args:
        df (pd.DataFrame): Dataframe containing image filenames and target data.
        filename (str): Column name for the image filenames.
        path (str): Path to the directory containing the images.
        train_split (float): Fraction of data to be used for training.
        val_split (float): Fraction of data to be used for validation.
        test_split (float): Fraction of data to be used for testing.
    
    Returns:
        tuple: A tuple containing:
            - df (pd.DataFrame): Updated dataframe with an additional img_path column.
            - train_df (pd.DataFrame): Dataframe for training.
            - val_df (pd.DataFrame): Dataframe for validation.
            - test_df (pd.DataFrame): Dataframe for testing.
    """
    # Ensure filenames are strings
    df[filename] = df[filename].astype(str)
    
    # Add the img_path column
    def construct_path(fn):
        if isinstance(fn, str):
            return os.path.abspath(os.path.join(path, fn if "." in fn else fn + ".jpg"))
        else:
            print(f"Non-string filename: {fn}")
            return None
    
    df['img_path'] = df[filename].apply(construct_path)
    
    # Split the dataframe
    train_df = df.sample(frac=train_split, random_state=SEED)
    val_df = df.drop(train_df.index).sample(frac=val_split / (val_split + test_split), random_state=SEED)
    test_df = df.drop(train_df.index).drop(val_df.index)
    
    return df, train_df, val_df, test_df


# Validate image paths
def validate_image_paths(df):
    """
    Validates the paths in the img_path column to ensure all images can be opened.
    
    Args:
        df (pd.DataFrame): Dataframe with an img_path column containing paths to the images.
    
    Returns:
        None. Prints an error message if any image path is invalid.
    """
    for path in df['img_path']:
        if cv2.imread(path) is None:
            print(f'Invalid image file: {path}')

# Generate batches of tensor image data
def create_data_generators(classification, train_df, val_df, test_df, y_col, path_col="img_path", batch_size_train=32, batch_size_val=1):
    """
    Generates batches of tensor image data for training, validation, and testing.
    
    Args:
        classification (bool): Whether the task is classification (True) or regression (False).
        train_df (pd.DataFrame): Dataframe for training data.
        val_df (pd.DataFrame): Dataframe for validation data.
        test_df (pd.DataFrame): Dataframe for testing data.
        y_col (str or list): Column names for the target data.
        path_col (str, optional): Column name for image paths. Default is "img_path".
        batch_size_train (int, optional): Batch size for training data. Default is 32.
        batch_size_val (int, optional): Batch size for validation data. Default is 1.
    
    Returns:
        tuple: A tuple containing:
            - train_generator (ImageDataGenerator): Generator for training data.
            - val_generator (ImageDataGenerator): Generator for validation data.
            - test_generator (ImageDataGenerator): Generator for testing data.
    """
    mode = "categorical" if classification else "raw"
    image_size = IMAGE_RESIZE
    datagen = preprocessing.image.ImageDataGenerator(rescale=1./255.)

    # Data generator for training data
    train_gen = datagen.flow_from_dataframe(train_df, x_col=path_col, y_col=y_col, batch_size=batch_size_train, class_mode=mode, shuffle=False, target_size=(image_size, image_size))
    
    # Data generator for validation data
    val_gen = datagen.flow_from_dataframe(val_df, x_col=path_col, y_col=y_col, batch_size=batch_size_val, class_mode=mode, shuffle=False, target_size=(image_size, image_size))
    
    # Data generator for testing data
    test_gen = datagen.flow_from_dataframe(test_df, x_col=path_col, y_col=y_col, batch_size=1, class_mode=mode, shuffle=False, target_size=(image_size, image_size))
    
    return train_gen, val_gen, test_gen

# Create and compile a transfer learning model
def create_model(base_model_name, learning_rate, optimizer_name, dense_layer_size, dense_units, final_activation, loss, metrics, dropout=0.2):
    """
    Creates and compiles a transfer learning model using a specified base model.
    
    Args:
        base_model_name (str): Name of the pre-trained base model.
        learning_rate (float): Learning rate for the optimizer.
        optimizer_name (str): Name of the optimizer ("Adam", "SGD", "RMSprop").
        dense_layer_size (int): Number of units in the dense layer added between the base model and the classifier added on top of it. Set to 0 to exclude.
        dense_units (int): Number of units in the final dense layer (top-level classifier).
        final_activation (str): Activation function for the final layer. Equals to the number of classes for classification. Equals to the number of outcome variables for regression.
        loss (str): Loss function for model training.
        metrics (list): List of metrics for model evaluation.
        dropout (float, optional): Dropout rate. Default is 0.2. Set to 0 to exclude.
    
    Returns:
        tuple: A tuple containing:
            - base_model (Model): The base model with pre-trained weights.
            - model (Model): The compiled model.
    """
    
    # Select optimizer
    opt = getattr(optimizers, optimizer_name)(learning_rate=learning_rate)
    
    # Load the pre-trained base model
    base_model_cls = getattr(applications, base_model_name)
    base_model = base_model_cls(include_top=False, pooling=POOLING_AVERAGE, weights='imagenet', input_shape=(IMAGE_RESIZE, IMAGE_RESIZE, 3))
    base_model.trainable = False

    # Model architecture
    inputs = layers.Input(shape=(IMAGE_RESIZE, IMAGE_RESIZE, 3))
    x = base_model(inputs, training=False)
    x = layers.Flatten()(x)
    # Add dense layer if specified
    if dense_layer_size > 0:
        x = layers.Dense(dense_layer_size, activation='relu')(x)
    # Add dropout layer if specified
    if dropout > 0:
        x = layers.Dropout(dropout)(x)
    outputs = layers.Dense(dense_units, activation=final_activation, name="prediction")(x)
    
    # Compile the model
    model = models.Model(inputs, outputs)
    model.compile(optimizer=opt, loss=loss, metrics=metrics)
    model.summary()
        
    return base_model, model

# Train the model with optional fine-tuning
def train_model(classification, fine_tuning, model, model_path, train_gen, val_gen, initial_epoch=0):
    """
    Trains a compiled model, optionally fine-tuning it.
    
    Args:
        classification (bool): Whether the task is classification (True) or regression (False).
        fine_tuning (bool): Whether to fine-tune the model (True) or not (False).
        model (Model): The compiled model.
        model_path (str): Path to save the trained model.
        train_gen (ImageDataGenerator): Generator for training data.
        val_gen (ImageDataGenerator): Generator for validation data.
        initial_epoch (int, optional): Starting epoch for training. Default is 0.
    
    Returns:
        History: Training history.
    """
    # Callbacks
    early_stopping = callbacks.EarlyStopping(monitor='val_loss', mode='min', patience=EARLY_STOP_PATIENCE, restore_best_weights=True, verbose=1)
    model_checkpoint = callbacks.ModelCheckpoint(filepath=model_path, monitor='val_loss', save_best_only=True, mode='min', verbose=1)
    
    # Compute class weights if classification
    class_weights = compute_class_weight("balanced", classes=np.unique(train_gen.classes), y=train_gen.classes) if classification else None
    
    # Train the model
    if class_weights is not None:
        class_weight_dict = dict(enumerate(class_weights))
        history = model.fit(train_gen, epochs=NUM_EPOCHS + (NUM_FINE_TUNE_EPOCHS if fine_tuning else 0), validation_data=val_gen, initial_epoch=initial_epoch, class_weight=class_weight_dict, callbacks=[early_stopping, model_checkpoint])
    else:
        history = model.fit(train_gen, epochs=NUM_EPOCHS + (NUM_FINE_TUNE_EPOCHS if fine_tuning else 0), validation_data=val_gen, initial_epoch=initial_epoch, callbacks=[early_stopping, model_checkpoint])
    
    model.load_weights(model_path)
    print(f'Classification: {classification}, Fine tuning: {fine_tuning}, Batch size: {train_gen.batch_size}, Optimizer: {model.optimizer}, Learning rate: {model.optimizer.learning_rate.numpy()}, Dropout rate: {model.layers[-2].rate if isinstance(model.layers[-2], layers.Dropout) else 0.0}, Dense layer size: {model.layers[-2].units if isinstance(model.layers[-2], layers.Dense) else 0}')
    return history

def unfreeze_model(base_model, model, optimizer_name, lr, loss):
    """
    Unfreezes layers in the base model for fine-tuning and recompiles the model.
    
    Args:
        base_model (Model): The pre-trained base model.
        model (Model): The compiled model that has been trained by transfer learning. 
        optimizer_name (str): Name of the optimizer ("Adam", "SGD", "RMSprop").
        lr (float): Learning rate for fine-tuning.
        loss (str): Loss function for model training.
    
    Returns:
        Model: The recompiled model with unfrozen layers.
    """
    base_model.trainable = True
    for layer in base_model.layers[:100]:
        layer.trainable = False
    opt = getattr(optimizers, optimizer_name)(learning_rate=lr)
    model.compile(optimizer=opt, loss=loss, metrics=['accuracy', 'mae'])


# Predict the testing dataset
def predict_new(classification, data_gen, df, model, model_name, y_cols, result_path):
    """
    Make predictions on certain target dataset
    
    Args:
        classification (bool): Whether the task is classification (True) or regression (False).
        data_gen (ImageDataGenerator): Generator for the target data.
        df (pd.DataFrame): Dataframe for the target data.
        model (Model): The trained model.
        model_name (str): Name of the model.
        y_cols (str or list): Column names for the target data.
        result_path (str): Path to save the results.
    
    Returns:
        np.ndarray: Predictions on the target data.
    """
    # Save evaluation metrics
    data_gen.reset()
    evals = model.evaluate(data_gen)
    evals_df = pd.DataFrame(evals)
    evals_df.index = ['loss', 'metric']
    evals_df.to_csv(os.path.join(result_path, f"test_metrics_{model_name}.csv"))
    
    """
    metrics_df = pd.DataFrame({'Validation Loss': [evals[0]]})
    if classification:
        metrics_df['Validation Accuracy'] = [evals[1]]
    else:
        metrics_df['Validation MAE'] = [evals[1]]
    metrics_df.to_csv(os.path.join(result_path, f"eval_{model_name}.csv"), index=False)
    """
    
    # Save predictions
    data_gen.reset()
    scores = model.predict(data_gen)
    result_df = pd.DataFrame(scores)
    if classification:
        result_df.columns = list(data_gen.class_indices.keys())
        result_df.insert(0, "label.pred", np.argmax(scores, axis=1))
    else:
        result_df.columns = y_cols
    
    fn = [os.path.basename(obj).replace(".jpg", "") for obj in data_gen.filenames]
    result_df.insert(0, "X", fn)
    result_df.to_csv(os.path.join(result_path, f"test_pred_{model_name}.csv"), index=False)
    
    # Save classification summary if classification
    if classification:
        report = skmetrics.classification_report(data_gen.classes, np.argmax(scores, axis=1), target_names=data_gen.class_indices)
        confusion = skmetrics.confusion_matrix(data_gen.classes, np.argmax(scores, axis=1), labels=[0, 1])
        with open(os.path.join(result_path, f"report_{model_name}.txt"), 'w') as f:
            f.write(f"Classification report:\n{report}\nConfusion matrix:\n{np.array2string(confusion, separator=', ')}")
    return scores


# Plot learning curves
def plot_learning_curves(classification, history, title, fine_tuning):
    """
    Plots learning curves for training and validation data and saves metrics.
    
    Args:
        classification (bool): Whether the task is classification (True) or regression (False).
        history (History): Training history.
        title (str): Title for the plot.
        fine_tuning (bool): Whether the model has been fine-tuned (True) or not (False).
    
    Returns:
        None. Saves and displays the plots.
    """
    
    history_df = pd.DataFrame(history.history)
    with open(os.path.join(EVAL_PATH, f'metrics_history_{title}.csv'), mode='w') as f:
        history_df.to_csv(f, index=False)
    
    metric = list(history.history.keys())[1]
    val_metric = list(history.history.keys())[3]
    
    plt.figure(figsize=(15, 8))
    plt.subplot(221)
    plt.plot(history.history[metric])
    plt.plot(history.history[val_metric])
    if fine_tuning:
        plt.axvline(x=NUM_EPOCHS, color='r', linestyle='--', label='Fine Tuning Start')
    plt.title(f'Model {metric.capitalize()}')
    plt.ylabel(metric.capitalize())
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Validation'])
    
    plt.subplot(222)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    if fine_tuning:
        plt.axvline(x=NUM_EPOCHS, color='r', linestyle='--', label='Fine Tuning Start')
    plt.title('Model Loss')
    plt.ylabel('MSE')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Validation'])
    
    plt.suptitle(title, fontsize=16)
    plt.savefig(os.path.join(PLOT_PATH, f"{title.replace(' ', '_')}.png"))
    plt.show()

# Plot ROC curve
def plot_roc(y_true, y_pred, target_names, title):
    """
    Plots the ROC curve for classification task.
    
    Args:
        y_true (np.ndarray): True labels.
        y_pred (np.ndarray): Predicted labels.
        target_names (list): List of target class names.
        title (str): Title for the plot.
    
    Returns:
        None. Saves and displays the ROC plot.
    """
    n_classes = len(target_names)
    y_true = utils.to_categorical(y_true, num_classes=n_classes)
    
    fpr, tpr, roc_auc = dict(), dict(), dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = skmetrics.roc_curve(y_true[:, i], y_pred[:, i])
        roc_auc[i] = skmetrics.auc(fpr[i], tpr[i])
    
    # Compute micro-average ROC curve and ROC area
    fpr['micro'], tpr['micro'], _ = skmetrics.roc_curve(y_true.ravel(), y_pred.ravel())
    roc_auc['micro'] = skmetrics.auc(fpr['micro'], tpr['micro'])
    
    # Compute macro-average ROC curve and ROC area
    # Aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
    # Interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
    # Average it and compute AUC
    mean_tpr /= n_classes
    
    fpr['macro'] = all_fpr
    tpr['macro'] = mean_tpr
    roc_auc['macro'] = skmetrics.auc(fpr['macro'], tpr['macro'])

    plt.figure(figsize=(10, 10))
    plt.plot(fpr['micro'], tpr['micro'], label=f'micro-average ROC curve (area = {roc_auc["micro"]:.2f})', color='deeppink', linestyle=':', linewidth=4)
    plt.plot(fpr['macro'], tpr['macro'], label=f'macro-average ROC curve (area = {roc_auc["macro"]:.2f})', color='navy', linestyle=':', linewidth=4)
    
    colors = itertools.cycle(['aqua', 'darkorange', 'cornflowerblue'])
    for i, color in zip(range(n_classes), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=2, label=f'ROC curve of class {target_names[i]} (area = {roc_auc[i]:.2f})')
    
    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=13)
    plt.ylabel('True Positive Rate', fontsize=13)
    plt.title(title, fontsize=15, pad=10)
    plt.legend(loc="lower right", fontsize=11)
    plt.savefig(os.path.join(ROC_PATH, f"{title.replace(' ', '_')}.png"))
    plt.show()
    print("ROC plot saved in " + os.path.join(ROC_PATH, f"{title.replace(' ', '_')}.png") + ".")
    
# Function to run hyperparameters tuning and model training
def run_hyperparams(base_model, classification, fine_tuning, train_df, val_df, test_df, batch_size, lr, ft_lr, opt, dropout, dense_size, y_cols, name, img_path_col='img_path'):
    """
    Runs hyperparameter tuning and trains the model.
    
    Args:
        base_model (str): Name of the pre-trained base model.
        classification (bool): Whether the task is classification (True) or regression (False).
        fine_tuning (bool): Whether to fine-tune the model (True) or not (False).
        train_df (pd.DataFrame): Dataframe for training data.
        val_df (pd.DataFrame): Dataframe for validation data.
        test_df (pd.DataFrame): Dataframe for test data.
        batch_size (int): Batch size for training.
        lr (float): Learning rate for training.
        ft_lr (float): Learning rate for fine-tuning.
        opt (str): Optimizer for training.
        dropout (float): Dropout rate.
        dense_size (int): Number of units in the dense layer.
        y_cols (str or list): Column names for the target data.
        name (str): Name of the model.
        img_path_col (str, optional): Column name for image paths. Default is 'img_path'.
    
    Returns:
        Model: The trained model.
    """
    train_gen, val_gen, test_gen = create_data_generators(classification, train_df, val_df, test_df, y_cols, img_path_col, batch_size, BATCH_SIZE_VALIDATION)
    
    dense_unit = len(train_gen.class_indices) if classification else len(y_cols)
    final_activation = "softmax" if classification or dense_unit > 1 else "linear"
    loss = "categorical_crossentropy" if classification else "mean_squared_error"
    metrics = ['accuracy'] if classification else ['mae']
    
    base_model, model = create_model(base_model, lr, opt, dense_size, dense_unit, final_activation, loss, metrics, dropout)
    print(f'Classification: {classification}, Fine tuning: {fine_tuning}, Batch size: {batch_size}, Optimizer: {opt}, Learning rate: {lr}, Dropout rate: {dropout}, Dense layer size: {dense_size}')
    print(f'Images are obtained from: {IMAGE_PATH}')
    
    model_path = os.path.join(MODEL_PATH, f"{name}.hdf5")
    
    history = train_model(classification, fine_tuning, model, model_path, train_gen, val_gen)
    plot_learning_curves(classification, history, name, fine_tuning)
    
    # scores = predict_new(classification, val_gen, val_df, model, name, y_cols, EVAL_PATH)
    test_scores = predict_new(classification, test_gen, test_df, model, name, y_cols, EVAL_PATH)
    
    if classification:
        plot_roc(test_gen.classes, test_scores, list(test_gen.class_indices.keys()), f'ROC curve_test_{name}')
        # plot_roc(val_gen.classes, scores, list(val_gen.class_indices.keys()), f'ROC curve_val_{name}')
    
    if fine_tuning:
        unfreeze_model(base_model, model, opt, ft_lr, loss)
        fine_name = f"{name}_ft_{ft_lr}"
        fine_model_path = os.path.join(MODEL_PATH, f"{fine_name}.hdf5")
        fine_history = train_model(classification, fine_tuning, model, fine_model_path, train_gen, val_gen, initial_epoch=history.epoch[-1])
        plot_learning_curves(classification, fine_history, fine_name, fine_tuning)


def run(project, task, data_file, outcome_list, patch_id, image_path, result_path, base_model, image, optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, test_split, val_split, seed, image_resize, num_epoch, patience):
    """
    Main function to run the training process with specified parameters.
    
    Args:
        project (str): Name of the project.
        task (str): Type of task ("classification" or "regression").
        data_file (str): Path to the CSV file containing the data.
        outcome (str or list): Column name(s) for the target variable(s).
        patch_id (str): Column name for the image IDs.
        image_path (str): Path to the directory containing the images.
        result_path (str): Path to the directory to save results.
        base_model (str): Name of the base model for transfer learning.
        image (str): Type of image normalization ("o", "m", "v").
        optimizer (str): Optimizer for training.
        batch_size (int): Batch size for training.
        learning_rate (float): Learning rate for training.
        dropout_rate (float): Dropout rate for the dropout layer.
        dense_layer_size (int): Number of units in the dense layer.
        test_split (float): Ratio for the test dataset split.
        val_split (float): Ratio for the validation dataset split.
    
    Returns:
        None.
    """
    global SEED, IMAGE_RESIZE, NUM_EPOCHS, EARLY_STOP_PATIENCE, IMAGE_PATH
    SEED = seed
    IMAGE_RESIZE = image_resize
    NUM_EPOCHS = num_epoch
    EARLY_STOP_PATIENCE = patience
    IMAGE_PATH = image_path

    model_name = f"{base_model}_{image}_{optimizer}_{batch_size}_{learning_rate}_{dropout_rate}_{dense_layer_size}"    

    
    if len(outcome_list) == 1 and task == "classification":
        outcome_list = outcome_list[0]
        
    setup_paths(result_path, project, outcome_list, task)
    
    if os.path.isfile(os.path.join(EVAL_PATH, f"test_pred_{model_name}.csv")):
        sys.exit()
    
    df = pd.read_csv(data_file)
    if patch_id not in df.columns:
        if 'Unnamed: 0' in df.columns:
            df.rename(columns={'Unnamed: 0': patch_id}, inplace=True)
        else:
            raise ValueError(f"Column {patch_id} does not exist in the dataframe. Available columns are: {df.columns.tolist()}")
    # Ensure 'Unnamed: 0' is not in columns
    if 'Unnamed: 0' in df.columns:
        df.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Ensure patch_id column exists after renaming
    if patch_id not in df.columns:
        raise ValueError(f"Column {patch_id} does not exist in the dataframe after renaming. Available columns are: {df.columns.tolist()}")


    if task == "classification":
        df[outcome_list] = df[outcome_list].apply(str)
        is_classification = True
    else:
        df[outcome_list] = df[outcome_list].apply(pd.to_numeric)
        is_classification = False
    
    df, train_df, val_df, test_df = prepare_dataframe(df, patch_id, image_path, 1 - test_split - val_split, val_split, test_split)
    outcome_str = '_'.join(outcome_list) if isinstance(outcome_list, list) else outcome_list
    train_df.to_csv(os.path.join(result_path, f'training_{project}_{outcome_str}.csv'))
    val_df.to_csv(os.path.join(result_path, f'validation_{project}_{outcome_str}.csv'))
    test_df.to_csv(os.path.join(result_path, f'testing_{project}_{outcome_str}.csv'))
    
    run_hyperparams(base_model, is_classification, False, train_df, val_df, test_df, batch_size, learning_rate, 0.000001, optimizer, dropout_rate, dense_layer_size, outcome_list, model_name)
    
   
def parse_args():
    parser = argparse.ArgumentParser(description='Take inputs')
    parser.add_argument('--project', default="unfilt", type=str, help='a string of project name')
    parser.add_argument('--task', default="classification", type=str, help='a string of task type (classification, regression)')
    parser.add_argument('--data_file', default="BRCA_clusters_for_STpath.csv", type=str, help='a string of the name of the csv file containing the predictor and response variables')
    parser.add_argument('--result_path', default="../output/BRCA/Classification/", type=str, help='a string of the directory saving all results of the project')
    parser.add_argument('--image_path', default="../../st2image_data/BRCA/output/patch_jpg_2/", type=str, help='a string of the path to the image patches')
    parser.add_argument('-n', '--outcome_list', default=[], nargs='+', help='a list of column names for the response variable')
    parser.add_argument('--patch_id', default="X", type=str, help='a string of the column name of the ID of image patches')
    parser.add_argument('--base_model', default="ResNet50", type=str, help='a string of base model used for transfer learning (efficientnet.EfficientNetB0 to efficientnet.EfficientNetB7, VGG16, ResNet101, MobileNet, MobileNetV2)')
    parser.add_argument('--image', default="o", type=str, help='a string of stain normalization type of input patch (o = original, m = macenko, v = vahadane)')
    parser.add_argument('--optimizer', default="Adam", type=str, help='a string of optimizer')
    parser.add_argument('--batch_size', default=128, type=int, help='an integer of the batch size for training and validation')
    parser.add_argument('--learning_rate', default=0.0001, type=float, help='a number of learning rate')
    parser.add_argument('--dropout_rate', default=0.0, type=float, help='a number of dropout rate')
    parser.add_argument('--dense_layer_size', default=256, type=int, help='an integer of intermediate dense layer')
    parser.add_argument('--test_split', default=0.15, type=float, help='a number of the split ratio of testing dataset')
    parser.add_argument('--val_split', default=0.15, type=float, help='a number of the split ratio of validation dataset')
    parser.add_argument('--seed', default=17, type=int, help='an integer of random seed')
    parser.add_argument('--image_resize', default=224, type=int, help='an integer of the size of image input to the model')
    parser.add_argument('--num_epoch', default=500, type=int, help='an integer of the maximum number of epochs')
    parser.add_argument('--patience', default=20, type=int, help='an integer of the early stopping patience')
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()
    print(args)
    run(**vars(args))
    
    