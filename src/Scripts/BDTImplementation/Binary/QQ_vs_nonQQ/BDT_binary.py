import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import numpy as np
from sklearn.impute import SimpleImputer
import os
import argparse
import logging
import xgboost as xgb
import time

# Configure logging
def setup_logging(script_name, output_dir):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create handlers
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(os.path.join(output_dir, f"{script_name}.log"))

    # Create formatters and add them to handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

def main():
    parser = argparse.ArgumentParser(description="Train binary BDT model on parquet files.")
    parser.add_argument(
        '-t', '--timestamp',
        type=str,
        required=True,
        help="Specify the timestamp used in the parquet files."
    )
    args = parser.parse_args()

    timestamp = args.timestamp
    output_dir = f'outputs/{timestamp}'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    input_dir = '/home/mukund/Projects/updatedCoffea/Coffea_Analysis/src/Scripts/BDTImplementation/inputs'

    setup_logging(os.path.splitext(os.path.basename(__file__))[0], output_dir)

    era = 'UL2016preVFP'

    parquet_files = [f for f in os.listdir(input_dir) if f.endswith('.parquet') and f.startswith('events_')]
    logging.info(f"Reading {len(parquet_files)} parquet files...")

    df = pd.concat([pd.read_parquet(os.path.join(input_dir, f)) for f in parquet_files], ignore_index=True)

    df['y_binary'] = df['y'].apply(lambda x: 0 if x == 0 else 1)

    logging.info(f"Number of events in class 0: {len(df[df['y_binary'] == 0])}")
    logging.info(f"Number of events in class 1: {len(df[df['y_binary'] == 1])}")

    # Extract 500000 events from each class
    df_0_sampled = df[df['y_binary'] == 0].sample(n=500000, random_state=42, replace=True)
    df_1_sampled = df[df['y_binary'] == 1].sample(n=500000, random_state=42, replace=True)

    # Combine the sampled data
    df_sampled = pd.concat([df_0_sampled, df_1_sampled], ignore_index=True)

    # Split the data into 30% for testing and 70% for training
    df_train, df_test = train_test_split(df_sampled, test_size=0.3, random_state=42, stratify=df_sampled['y_binary'])

    logging.info(f"Number of events in training data for class 0: {len(df_train[df_train['y_binary'] == 0])}")
    logging.info(f"Number of events in training data for class 1: {len(df_train[df_train['y_binary'] == 1])}")

    # features = ['ttbarpz', 'ttbar_mass', 'Jet_HT', 'nJet', 'pT_sum', 'FW1', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'S', 'P', 'A', 'p2in', 'p2out']
    features = ['FW1', 'nJet', 'pT_sum', 'P', 'A', 'p2in', 'Syy', 'Sxy']
    logging.info(f"Features used for training: {features}")

    # Select features and target variable for training
    X_train_binary = df_train[features]
    y_train_binary = df_train['y_binary']

    # Select features and target variable for testing
    X_test_binary = df_test[features]
    y_test_binary = df_test['y_binary']

    # Check for NaNs in the original data
    logging.info("Checking for NaNs in the original data...")
    logging.info(f"NaNs in X_train_binary: {np.isnan(X_train_binary).sum()}")
    logging.info(f"NaNs in X_test_binary: {np.isnan(X_test_binary).sum()}")

    # Impute missing values
    imputer = SimpleImputer(strategy='mean')
    X_train_binary = imputer.fit_transform(X_train_binary)
    X_test_binary = imputer.transform(X_test_binary)

    # Check for remaining NaNs after imputation
    logging.info("Checking for NaNs after imputation...")
    logging.info(f"NaNs in X_train_binary: {np.isnan(X_train_binary).sum()}")
    logging.info(f"NaNs in X_test_binary: {np.isnan(X_test_binary).sum()}")

    # Set up the XGBoost model
    xgb_binary = xgb.XGBClassifier(tree_method='gpu_hist', random_state=42)

    # Measure the time to train a single model
    start_time = time.time()
    xgb_binary.fit(X_train_binary, y_train_binary)
    end_time = time.time()

    single_model_time = end_time - start_time
    logging.info(f"Single model training time: {single_model_time} seconds")

    # exit()

    # Define parameter grid for GridSearchCV
    param_grid = {
        'n_estimators': [50, 100, 200],
        'learning_rate': [0.01, 0.1, 0.2],
        'max_depth': [3, 4, 5],
        'min_child_weight': [1, 2, 4],
        'subsample': [0.8, 0.9, 1.0],
        'colsample_bytree': [0.8, 0.9, 1.0]
    }

    # Perform GridSearchCV
    # optimizing_para = 'roc_auc'
    optimizing_para = 'accuracy'
    logging.info(f"Optimizing parameter: {optimizing_para}")
    grid_search = GridSearchCV(estimator=xgb_binary, param_grid=param_grid, cv=3, scoring= optimizing_para, n_jobs=-1)
    logging.info("Performing GridSearchCV...")
    grid_search.fit(X_train_binary, y_train_binary)

    # Get the best parameters and train the final model
    best_params = grid_search.best_params_
    logging.info(f"Best parameters found: {best_params}")

    xgb_binary_best = xgb.XGBClassifier(**best_params, tree_method='gpu_hist', random_state=42)
    xgb_binary_best.fit(X_train_binary, y_train_binary)

    # Predict probabilities for the testing dataset
    y_pred_proba_binary = xgb_binary_best.predict_proba(X_test_binary)[:, 1]

    # Predict classes for the testing dataset
    y_pred_binary = xgb_binary_best.predict(X_test_binary)

    # Calculate accuracy
    accuracy = accuracy_score(y_test_binary, y_pred_binary)
    logging.info(f"Accuracy: {accuracy}")

    # Calculate ROC AUC score
    roc_auc = roc_auc_score(y_test_binary, y_pred_proba_binary)
    logging.info(f"ROC AUC Score: {roc_auc}")

    # Plot ROC curve
    fpr, tpr, _ = roc_curve(y_test_binary, y_pred_proba_binary)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.savefig(f'{output_dir}/roc_curve_binary.png')
    plt.close()
    logging.info(f"ROC curve saved to {output_dir}/roc_curve_binary.png")


    # Compute ROC curve for training and test sets
    fpr_train, tpr_train, _ = roc_curve(y_train_binary, xgb_binary_best.predict_proba(X_train_binary)[:, 1])
    fpr_test, tpr_test, _ = roc_curve(y_test_binary, y_pred_proba_binary)

    # Plot both training and test ROC curves
    plt.figure()
    plt.plot(fpr_train, tpr_train, color='blue', lw=2, label=f'Training ROC (area = {roc_auc_score(y_train_binary, xgb_binary_best.predict_proba(X_train_binary)[:, 1]):.2f})')
    plt.plot(fpr_test, tpr_test, color='darkorange', lw=2, label=f'Test ROC (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')  # Random guess line
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Training vs Test ROC Curve')
    plt.legend(loc="lower right")

    # Save the plot
    plt.savefig(f'{output_dir}/roc_curve_train_vs_test.png')
    logging.info(f"Training vs Test ROC curve saved to {output_dir}/roc_curve_train_vs_test.png")


    # Create a new DataFrame to store the actual y values and the output scores
    df_scores_binary = pd.DataFrame({
        'y_actual': y_test_binary,
        'score': y_pred_proba_binary
    })

    # Save the DataFrame to a CSV file
    df_scores_binary.to_csv(f'{output_dir}/scores_binary.csv', index=False)
    logging.info(f"Scores saved to {output_dir}/scores_binary.csv")

    # Optional: Save the model for later use (if needed)
    joblib.dump(xgb_binary_best, f'{output_dir}/xgb_model_binary.pkl')
    logging.info("Model saved")
    logging.info("="*100)

if __name__ == '__main__':
    main()
