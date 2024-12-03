import pandas as pd
import numpy as np
import logging
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

# Set up logging to file and console
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
file_handler = logging.FileHandler('output.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(file_handler)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console_handler)

# Load your data
data_0 = pd.read_parquet('inputs/events_UL2016preVFP_chunk0_AL.parquet')
data_1 = pd.read_parquet('inputs/events_UL2016preVFP_chunk1_AL.parquet')

# Combine the data
data = pd.concat([data_0, data_1], ignore_index=True)

# Select only the specified observables
all_observables = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'Jet_HT', 'S', 'P', 'A', 'p2in', 'p2out', 'AL', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']


y = data['y'].map({-1: 1, 0: 0, 1: 1})

# Dictionary to store last trained epoch info
last_trained_epoch_info = {}


def train_and_evaluate(observables, iteration):
    logger.info(f"Training with observables: {observables}")

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
    y_balanced = to_categorical(y_balanced, num_classes=2)

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
    model.add(Dense(2, activation='softmax'))
    # Compile the model
    model.compile(loss='categorical_crossentropy', optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), metrics=['accuracy'])

    # Early stopping
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

    # Train the model
    epochs = 50
    batchsize = 32
    validation_split = 0.2
    history = model.fit(X_train, y_train, epochs= epochs, batch_size= batchsize, validation_split= validation_split, callbacks=[early_stopping])

    # Evaluate the model
    y_pred = model.predict(X_test)
    y_pred_classes = np.argmax(y_pred, axis=1)
    y_test_classes = np.argmax(y_test, axis=1)

    accuracy =  accuracy_score(y_test_classes, y_pred_classes)

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
    last_epoch_info = {}
    for key in history.history:
        last_epoch_info[key] = history.history[key][-1]
        logger.info(f"{key}: {history.history[key][-1]}")
    last_trained_epoch_info[f'model{iteration}'] = last_epoch_info

    # Print classification report
    logger.info("Classification report:")
    logger.info(classification_report(y_test_classes, y_pred_classes))

    return accuracy


# Initial training with all observables
current_observables = all_observables.copy()
iteration = 0
best_accuracy = train_and_evaluate(current_observables, iteration)

# Iteratively remove the least effective feature
while len(current_observables) > 1:
    accuracies = {}
    for feature in current_observables:
        temp_observables = current_observables.copy()
        temp_observables.remove(feature)
        accuracy = train_and_evaluate(temp_observables, f"{iteration}_{feature}")
        accuracies[feature] = accuracy

    # Find the feature whose removal causes the smallest drop in accuracy
    least_effective_feature = min(accuracies, key=lambda k: accuracies[k])
    logger.info("="*30)
    logger.info(f"At the end of iteration {iteration}, the least effective feature is {least_effective_feature} with an accuracy of {accuracies[least_effective_feature]}")
    
    new_accuracy = accuracies[least_effective_feature]

    # Check if the new accuracy is significantly lower
    if new_accuracy < best_accuracy - 0.05:  # Adjust the threshold as needed
        break

    # Update the observables and best accuracy
    current_observables.remove(least_effective_feature)
    logger.info(f"New accuracy after removing {least_effective_feature}: {new_accuracy}")
    logger.info(f"Current observables after iteration {iteration}: {current_observables}")
    logger.info("="*30)
    # best_accuracy = new_accuracy
    iteration += 1

# Save the last trained epoch info to a file
with open('last_trained_epoch_info.json', 'w') as f:
    import json
    json.dump(last_trained_epoch_info, f, indent=4)