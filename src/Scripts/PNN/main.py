import numpy as np
import os
import pandas as pd
import tensorflow as tf
import argparse
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import shap

era = 'UL2016preVFP'
parser = argparse.ArgumentParser(description="Train BDT models on parquet files.")
parser.add_argument(
    '-t', '--timestamp',
    type=str,
    required=True,
    help="Specify the timestamp used in the parquet files."
)
args = parser.parse_args()

timestamp = args.timestamp
outputDir = f'outputs/{timestamp}'
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

parquet_files = [f for f in os.listdir(outputDir) if f.endswith('.parquet') and f.startswith('events_')]
df = pd.concat([pd.read_parquet(os.path.join(outputDir, f)) for f in parquet_files], ignore_index=True)


df_qq = df[df['y'] == 0]
df_gg = df[df['y'] == 1]
df_qg = df[df['y'] == 2]

# Extract 500000 events from each class
df_qq_sampled = df_qq.sample(n=500000, random_state=42, replace=True)
df_gg_sampled = df_gg.sample(n=500000, random_state=42, replace=True)
df_qg_sampled = df_qg.sample(n=500000, random_state=42, replace=True)


df = pd.concat([df_qq_sampled, df_gg_sampled, df_qg_sampled], ignore_index=True)

df_train, df_test = train_test_split(df, test_size=0.3, random_state=42, stratify=df['y'])

obss = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'P', 'A', 'p2in', 'Syy', 'Sxy']

# Select features and target variable for training
X_train = df_train[obss]
y_train = df_train['y']

# Select features and target variable for testing
X_test = df_test[obss]
y_test = df_test['y']


# 2. Preprocess Data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_train)

# Convert labels to one-hot encoding
n_classes = len(np.unique(y_train))
y_onehot = tf.keras.utils.to_categorical(y_train, num_classes=n_classes)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_onehot, test_size=0.2, random_state=42)

# 3. Define the Physics-Aware Neural Network Model
def physics_constraint_loss(y_true, y_pred):
    
    energy_pred = tf.reduce_sum(y_pred, axis=1)
    energy_true = tf.reduce_sum(y_true, axis=1)
    constraint_violation = tf.square(energy_pred - energy_true)
    return tf.reduce_mean(constraint_violation)


def custom_loss(y_true, y_pred):
    classification_loss = tf.keras.losses.categorical_crossentropy(y_true, y_pred)
    constraint_loss = physics_constraint_loss(y_true, y_pred)
    return classification_loss + 0.1 * constraint_loss  # Lambda controls constraint weight

# Build the model
model = Sequential([
    Dense(128, activation='relu', input_shape=(X_train.shape[1],)),
    BatchNormalization(),
    Dropout(0.5),
    Dense(64, activation='relu'),
    BatchNormalization(),
    Dropout(0.5),
    Dense(32, activation='relu'),
    Dense(n_classes, activation='softmax')
])

model.compile(optimizer=Adam(learning_rate=0.001),
              loss=custom_loss,
              metrics=['accuracy'])

# 4. Train the Model
callbacks = [
    EarlyStopping(patience=10, restore_best_weights=True),
    ReduceLROnPlateau(patience=5)
]

history = model.fit(
    X_train, y_train,
    validation_split=0.2,
    epochs=100,
    batch_size=128,
    callbacks=callbacks
)

# 5. Evaluate the Model
results = model.evaluate(X_test, y_test, verbose=0)
print(f"Test Loss: {results[0]:.4f}, Test Accuracy: {results[1]:.4f}")

# 6. Save the Model
model.save("physics_aware_model.h5")

# # 7. Interpretation (Optional)
# # Use SHAP or LIME for feature importance analysis

# explainer = shap.Explainer(model, X_train)
# shap_values = explainer(X_test)
# shap.summary_plot(shap_values, X_test)
