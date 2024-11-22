import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report
from sklearn.utils import class_weight
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, LeakyReLU
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2
from tensorflow.keras.utils import to_categorical
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data
data_0 = pd.read_parquet('inputs/events_UL2016preVFP_chunk0_AL.parquet')
data_1 = pd.read_parquet('inputs/events_UL2016preVFP_chunk1_AL.parquet')

# Combine the data
data = pd.concat([data_0, data_1], ignore_index=True)

# Select only the specified observables
observables = ['ttbarpz', 'FW1', 'AL', 'P', 'p2in']
X = data[observables]
y = data['y']

y = y.map({-1: 1, 0: 0, 1: 1})

# Combine X and y for undersampling
data_balanced = pd.concat([X, y], axis=1)

# Find the minimum class count
min_class_count = data_balanced['y'].value_counts().min()

# Undersample each class to the minimum class count
data_balanced = data_balanced.groupby('y').apply(lambda x: x.sample(min_class_count)).reset_index(drop=True)

for observable in observables:
    plt.figure(figsize=(10, 6))
    sns.histplot(data=data_balanced, x=observable, hue='y', kde=True, element='step')
    plt.title(f'Distribution of {observable} for different y values')
    plt.savefig(f'outputs/{observable}.png')
    plt.close()

# Separate features and target
X_balanced = data_balanced.drop('y', axis=1)
y_balanced = data_balanced['y']

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
batchsize = 64
validation_split = 0.2
history = model.fit(X_train, y_train, epochs= epochs, batch_size= batchsize, validation_split= validation_split, callbacks=[early_stopping])

# Evaluate the model
y_pred = model.predict(X_test)
y_pred_classes = np.argmax(y_pred, axis=1)
y_test_classes = np.argmax(y_test, axis=1)

print("List of observables being used:")
print(observables)

print("Class distribution after undersampling:")
print(Nsamples)

# Print model summary
print("Model structure:")
model.summary()

# Print training parameters
print("Training parameters:")
print(f"Number of epochs: {epochs}")
print(f"Batch size: {batchsize}")
print(f"Validation split: {validation_split}")

# Print training history
print("Last trained epoch:")
for key in history.history:
    print(key, history.history[key][-1])

# Print classification report
print("Classification report:")
print(classification_report(y_test_classes, y_pred_classes))