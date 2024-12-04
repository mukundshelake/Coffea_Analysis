import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import numpy as np
from sklearn.impute import SimpleImputer

timestamp = 'Nov11'
outputDir = f'outputs/{timestamp}'
era = 'UL2016preVFP'
parquetFile = f"events_{era}.parquet"

# Load the Parquet file
df = pd.read_parquet(f"{outputDir}/{parquetFile}")


# Sample a fraction of the DataFrame
# df = df.sample(frac=0.1, random_state=42)
print(df.head())

# Split the data into 30% for testing and 70% for training
df_train, df_test = train_test_split(df, test_size=0.3, random_state=42, stratify=df['y'])

# Select features and target variable for training
X_train = df_train[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_train = df_train['y']

# Select features and target variable for testing
X_test = df_test[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_test = df_test['y']

# Prepare data for class 0 vs class 1
df_train_0_vs_1 = df_train[df_train['y'].isin([0, 1])]
X_train_0_vs_1 = df_train_0_vs_1[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_train_0_vs_1 = df_train_0_vs_1['y']

df_test_0_vs_1 = df_test[df_test['y'].isin([0, 1])]
X_test_0_vs_1 = df_test_0_vs_1[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_test_0_vs_1 = df_test_0_vs_1['y']

# Prepare data for class 0 vs class 2
df_train_0_vs_2 = df_train[df_train['y'].isin([0, -1])]
X_train_0_vs_2 = df_train_0_vs_2[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_train_0_vs_2 = df_train_0_vs_2['y']

df_test_0_vs_2 = df_test[df_test['y'].isin([0, -1])]
X_test_0_vs_2 = df_test_0_vs_2[['ttbarpz', 'FW1', 'FW2', 'FW3', 'FW4', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'nJet', 'pT_sum', 'Jet_HT']]
y_test_0_vs_2 = df_test_0_vs_2['y']


# Check for NaNs in the original data
print("Checking for NaNs in the original data...")
print(f"NaNs in X_train_0_vs_1: {np.isnan(X_train_0_vs_1).sum()}")
print(f"NaNs in X_test_0_vs_1: {np.isnan(X_test_0_vs_1).sum()}")
print(f"NaNs in X_train_0_vs_2: {np.isnan(X_train_0_vs_2).sum()}")
print(f"NaNs in X_test_0_vs_2: {np.isnan(X_test_0_vs_2).sum()}")
print(f"NaNs in X_test: {np.isnan(X_test).sum()}")
# Impute missing values
imputer = SimpleImputer(strategy='mean')
X_train_0_vs_1 = imputer.fit_transform(X_train_0_vs_1)
X_test_0_vs_1 = imputer.transform(X_test_0_vs_1)
X_train_0_vs_2 = imputer.fit_transform(X_train_0_vs_2)
X_test_0_vs_2 = imputer.transform(X_test_0_vs_2)
X_test = imputer.transform(X_test)


# Check for remaining NaNs after imputation
print("Checking for NaNs after imputation...")
print(f"NaNs in X_train_0_vs_1: {np.isnan(X_train_0_vs_1).sum()}")
print(f"NaNs in X_test_0_vs_1: {np.isnan(X_test_0_vs_1).sum()}")
print(f"NaNs in X_train_0_vs_2: {np.isnan(X_train_0_vs_2).sum()}")
print(f"NaNs in X_test_0_vs_2: {np.isnan(X_test_0_vs_2).sum()}")
print(f"NaNs in X_test: {np.isnan(X_test).sum()}")
# Combine the testing data for prediction
X_test_combined = np.concatenate((X_test_0_vs_1, X_test_0_vs_2), axis=0)

# Check for NaNs in the combined testing data
print("Checking for NaNs in the combined testing data...")
print(f"NaNs in X_test_combined: {np.isnan(X_test_combined).sum()}")

# Set up the BDT models
bdt_0_vs_1 = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
bdt_0_vs_2 = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)

# Train the models
bdt_0_vs_1.fit(X_train_0_vs_1, y_train_0_vs_1)
bdt_0_vs_2.fit(X_train_0_vs_2, y_train_0_vs_2)

# Predict probabilities for the testing dataset
y_pred_proba_0_vs_1 = bdt_0_vs_1.predict_proba(X_test)[:, 1]
y_pred_proba_0_vs_2 = bdt_0_vs_2.predict_proba(X_test)[:, 1]

# Create a new DataFrame to store the actual y values and the two output scores
df_scores = pd.DataFrame({
    'y_actual': y_test,
    'score_0_vs_1': y_pred_proba_0_vs_1,
    'score_0_vs_2': y_pred_proba_0_vs_2
})

print(df_scores.head())

# Save the DataFrame to a CSV file
df_scores.to_csv(f'{outputDir}/scores.csv', index=False)

# Plot the 2D scatter plot
plt.figure(figsize=(10, 6))
scatter = plt.scatter(df_scores['score_0_vs_1'], df_scores['score_0_vs_2'], c=df_scores['y_actual'], cmap='viridis', alpha=0.6)
plt.colorbar(scatter, label='Actual Class')
plt.xlabel('Score for Class 0 vs Class 1')
plt.ylabel('Score for Class 0 vs Class 2')
plt.title('2D Scatter Plot of BDT Output Scores')
plt.savefig(f'{outputDir}/scatter_plot.png')
plt.show()

# Optional: Save the models for later use (if needed)
joblib.dump(bdt_0_vs_1, 'bdt_model_0_vs_1.pkl')
joblib.dump(bdt_0_vs_2, 'bdt_model_0_vs_2.pkl')


