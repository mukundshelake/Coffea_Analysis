import pandas as pd
import numpy as np
import logging
import json
import os
import argparse
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, accuracy_score
from sklearn.utils import class_weight
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, LeakyReLU
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2
from tensorflow.keras.utils import to_categorical
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

def setup_logging(output_dir):
    # Set up logging to file and console
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    file_handler = logging.FileHandler(os.path.join(output_dir, 'log_3classes.log'))
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(console_handler)
    return logger

def train_and_evaluate(observables, output_dir, logger):
    logger.info(f"Training with observables: {observables}")

    # Load your data
    data_0 = pd.read_parquet('inputs/events_UL2016preVFP_chunk0_AL.parquet')
    data_1 = pd.read_parquet('inputs/events_UL2016preVFP_chunk1_AL.parquet')

    # Combine the data
    data = pd.concat([data_0, data_1], ignore_index=True)

    # Select only the specified observables
    all_observables = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'Jet_HT', 'S', 'P', 'A', 'p2in', 'p2out', 'AL', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']

    # Map 'y' to three classes
    y = data['y'].map({-1: 2, 0: 0, 1: 1})

    # Select the specified observables
    X = data[list(observables)] 

    # Combine X and y for undersampling
    data_balanced = pd.concat([X, y], axis=1)

    # Find the minimum class count
    min_class_count = data_balanced['y'].value_counts().min()

    # Undersample each class to the minimum class count
    data_balanced = data_balanced.groupby('y').apply(lambda x: x.sample(min_class_count)).reset_index(drop=True)

    # Separate features and target
    X_balanced = data_balanced.drop('y', axis=1)
    y_balanced = data_balanced['y']

    # Get the class distribution after undersampling
    Nsamples = y_balanced.value_counts()

    # One-hot encode the target variable
    y_balanced = to_categorical(y_balanced, num_classes=3)

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_balanced, y_balanced, test_size=0.2, random_state=42)

    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Build the DNN model
    model = Sequential()
    model.add(Dense(128, input_shape=(X_train.shape[1],), activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(3, activation='softmax'))  # Update to 3 output nodes

    # Compile the model
    model.compile(loss='categorical_crossentropy', optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), metrics=['accuracy'])

    # Early stopping
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

    # Train the model
    epochs = 50
    batchsize = 32
    validation_split = 0.2
    history = model.fit(X_train, y_train, epochs=epochs, batch_size=batchsize, validation_split=validation_split, callbacks=[early_stopping])

    # Evaluate the model
    y_pred = model.predict(X_test)
    y_pred_classes = np.argmax(y_pred, axis=1)
    y_test_classes = np.argmax(y_test, axis=1)

    accuracy = accuracy_score(y_test_classes, y_pred_classes)

    logger.info("Class distribution after undersampling:")
    logger.info(Nsamples)

    # Print model summary
    logger.info("Model structure:")
    model.summary(print_fn=logger.info)

    # Print training parameters
    logger.info("Training parameters:")
    logger.info(f"Number of epochs: {epochs}")
    logger.info(f"Batch size: {batchsize}")
    logger.info(f"Validation split: {validation_split}")

    # Print training history
    logger.info("Last trained epoch:")
    for key in history.history:
        logger.info(f"{key}: {history.history[key][-1]}")

    # Print classification report
    logger.info("Classification report:")
    logger.info(classification_report(y_test_classes, y_pred_classes))

    # Save the outputs of the final layer's 3 nodes and the actual y values for the test sample
    output_data = {
        'y_test': y_test_classes.tolist(),
        'y_pred': y_pred.tolist()
    }
    with open(os.path.join(output_dir, 'output_data_3classes.json'), 'w') as f:
        json.dump(output_data, f, indent=4)

    return accuracy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train and evaluate DNN classifier.')
    parser.add_argument('-o', '--output', type=str, default='outputs', help='Output directory for the results.')
    parser.add_argument('-t', '--timestamp', type=str, default='timestamp', help='Timestamp to create a subfolder for the results.')

    args = parser.parse_args()

    output_dir = os.path.join(args.output, args.timestamp)
    os.makedirs(output_dir, exist_ok=True)

    logger = setup_logging(output_dir)

    # Initial training with all observables
    current_observables = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'Jet_HT', 'S', 'P', 'A', 'p2in', 'p2out', 'AL', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']
    best_accuracy = train_and_evaluate(current_observables, output_dir, logger)