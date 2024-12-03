import pandas as pd
import logging
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, LeakyReLU
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.utils import to_categorical
from scikeras.wrappers import KerasClassifier
import os

# Set up logging to file and console
script_name = os.path.splitext(os.path.basename(__file__))[0]
log_filename = f"{script_name}.log"
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
file_handler = logging.FileHandler(log_filename)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(file_handler)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console_handler)

logger.info("Starting the script...")

# Load your data
logger.info("Loading data...")
data_0 = pd.read_parquet('inputs/events_UL2016preVFP_chunk0_AL.parquet')
data_1 = pd.read_parquet('inputs/events_UL2016preVFP_chunk1_AL.parquet')

# Combine the data
logger.info("Combining data...")
data = pd.concat([data_0, data_1], ignore_index=True)

# Select only the specified observables
all_observables = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'P', 'A', 'p2in', 'Syy', 'Sxy']
logger.info(f"Selected observables: {all_observables}")

y = data['y'].map({-1: 1, 0: 0, 1: 1})

# Dictionary to store last trained epoch info
last_trained_epoch_info = {}

def create_model(input_shape, learning_rate=0.001, dropout_rate=0.5, neurons_layer1=128, neurons_layer2=64, neurons_layer3=32, activation='relu', optimizer='adam'):
    model = Sequential()
    model.add(Dense(neurons_layer1, input_shape=input_shape, activation=activation))
    model.add(Dropout(dropout_rate))
    model.add(Dense(neurons_layer2, activation=activation))
    model.add(Dropout(dropout_rate))
    model.add(Dense(neurons_layer3, activation=activation))
    model.add(Dense(2, activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    return model

def train_and_evaluate(observables):
    logger.info(f"Training with observables: {observables}")

    # Select the specified observables
    X = data[list(observables)] 

    # Combine X and y for undersampling
    data_balanced = pd.concat([X, y], axis=1)

    # Find the minimum class count
    min_class_count = data_balanced['y'].value_counts().min()
    logger.info(f"Minimum class count for undersampling: {min_class_count}")

    # Undersample each class to the minimum class count
    data_balanced = data_balanced.groupby('y').apply(lambda x: x.sample(min_class_count)).reset_index(drop=True)

    # Separate features and target
    X_balanced = data_balanced.drop('y', axis=1)
    y_balanced = data_balanced['y']

    # Get the class distribution after undersampling
    Nsamples = y_balanced.value_counts()
    logger.info(f"Class distribution after undersampling: {Nsamples}")

    # One-hot encode the target variable
    y_balanced = to_categorical(y_balanced, num_classes=2)

    # Split the data into training and testing sets
    logger.info("Splitting data into training and testing sets...")
    X_train, X_test, y_train, y_test = train_test_split(X_balanced, y_balanced, test_size=0.2, random_state=42)

    # Standardize the features
    logger.info("Standardizing features...")
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Define early stopping callback
    early_stopping = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)

    # Wrap your model using KerasClassifier
    model = KerasClassifier(
        model=create_model,
        input_shape=(X_train.shape[1],),
        verbose=0
    )

    # Define the parameter grid
    param_grid = {
        'batch_size': [10, 20, 40],
        'epochs': [10, 50, 100],
        'model__learning_rate': [0.001, 0.01, 0.1],
        'model__dropout_rate': [0.3, 0.5, 0.7],
        'model__neurons_layer1': [64, 128, 256],
        'model__neurons_layer2': [32, 64, 128],
        'model__neurons_layer3': [16, 32, 64],
        'model__activation': ['relu', 'tanh'],
        'model__optimizer': ['adam', 'rmsprop', 'sgd']
    }
    logger.info(f"Parameter grid: {param_grid}")

    # Perform Grid Search
    logger.info("Starting Grid Search...")
    grid = GridSearchCV(estimator=model, param_grid=param_grid, n_jobs=-1, cv=3)
    grid_result = grid.fit(X_train, y_train, callbacks=[early_stopping])

    # Print the best parameters and best score
    logger.info(f"Best: {grid_result.best_score_} using {grid_result.best_params_}")

    # Train the final model with the best parameters
    best_params = grid_result.best_params_
    logger.info(f"Training final model with best parameters: {best_params}")
    final_model = create_model(
        input_shape=(X_train.shape[1],),
        learning_rate=best_params['model__learning_rate'],
        dropout_rate=best_params['model__dropout_rate'],
        neurons_layer1=best_params['model__neurons_layer1'],
        neurons_layer2=best_params['model__neurons_layer2'],
        neurons_layer3=best_params['model__neurons_layer3'],
        activation=best_params['model__activation'],
        optimizer=best_params['model__optimizer']
    )
    final_model.fit(X_train, y_train, epochs=best_params['epochs'], batch_size=best_params['batch_size'], validation_split=0.2, callbacks=[early_stopping], verbose=1)

    # Evaluate the final model
    logger.info("Evaluating final model...")
    accuracy = final_model.evaluate(X_test, y_test, verbose=0)
    logger.info(f"Final model accuracy: {accuracy[1]}")

    return accuracy[1]

# Train and evaluate with all observables
logger.info("Starting training and evaluation with all observables...")
best_accuracy = train_and_evaluate(all_observables)
logger.info(f"Best accuracy: {best_accuracy}")

# Save the last trained epoch info to a file
logger.info("Saving last trained epoch info to file...")
with open('last_trained_epoch_info.json', 'w') as f:
    import json
    json.dump(last_trained_epoch_info, f, indent=4)

logger.info("Script finished.")
