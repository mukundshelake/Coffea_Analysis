import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import joblib
import numpy as np
from sklearn.impute import SimpleImputer
import os
import argparse
import logging

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

    setup_logging(os.path.splitext(os.path.basename(__file__))[0], outputDir)

    era = 'UL2016preVFP'

    parquet_files = [f for f in os.listdir(outputDir) if f.endswith('.parquet') and f.startswith('events_')]
    logging.info(f"Reading {len(parquet_files)} parquet files...")

    df = pd.concat([pd.read_parquet(os.path.join(outputDir, f)) for f in parquet_files], ignore_index=True)

    df_qq = df[df['y'] == 0]
    df_gg = df[df['y'] == 1]
    df_qg = df[df['y'] == 2]

    logging.info(f"Number of events in class 0: {len(df_qq)}")
    logging.info(f"Number of events in class 1: {len(df_gg)}")
    logging.info(f"Number of events in class 2: {len(df_qg)}")

    # Extract 500000 events from each class
    df_qq_sampled = df_qq.sample(n=500000, random_state=42, replace=True)
    df_gg_sampled = df_gg.sample(n=500000, random_state=42, replace=True)
    df_qg_sampled = df_qg.sample(n=500000, random_state=42, replace=True)

    # Combine the sampled data
    df = pd.concat([df_qq_sampled, df_gg_sampled, df_qg_sampled], ignore_index=True)

    # Split the data into 30% for testing and 70% for training
    df_train, df_test = train_test_split(df, test_size=0.3, random_state=42, stratify=df['y'])

    logging.info(f"Number of events in training data for qq: {len(df_train[df_train['y'] == 0])}")
    logging.info(f"Number of events in training data for gg: {len(df_train[df_train['y'] == 1])}")
    logging.info(f"Number of events in training data for qg: {len(df_train[df_train['y'] == 2])}")

    obss = ['ttbarpz', 'FW1', 'nJet', 'pT_sum', 'P', 'A', 'p2in', 'Syy', 'Sxy']

    # Select features and target variable for training
    X_train = df_train[obss]
    y_train = df_train['y']

    # Select features and target variable for testing
    X_test = df_test[obss]
    y_test = df_test['y']

    # Prepare data for class 0 vs class 1
    df_train_0_vs_1 = df_train[df_train['y'].isin([0, 1])]
    X_train_0_vs_1 = df_train_0_vs_1[obss]
    y_train_0_vs_1 = df_train_0_vs_1['y']

    df_test_0_vs_1 = df_test[df_test['y'].isin([0, 1])]
    X_test_0_vs_1 = df_test_0_vs_1[obss]
    y_test_0_vs_1 = df_test_0_vs_1['y']

    # Prepare data for class 0 vs class 2
    df_train_0_vs_2 = df_train[df_train['y'].isin([0, 2])]
    X_train_0_vs_2 = df_train_0_vs_2[obss]
    y_train_0_vs_2 = df_train_0_vs_2['y']

    df_test_0_vs_2 = df_test[df_test['y'].isin([0, 2])]
    X_test_0_vs_2 = df_test_0_vs_2[obss]
    y_test_0_vs_2 = df_test_0_vs_2['y']

    # Prepare data for class 1 vs class 2
    df_train_1_vs_2 = df_train[df_train['y'].isin([1, 2])]
    X_train_1_vs_2 = df_train_1_vs_2[obss]
    y_train_1_vs_2 = df_train_1_vs_2['y']

    df_test_1_vs_2 = df_test[df_test['y'].isin([1, 2])]
    X_test_1_vs_2 = df_test_1_vs_2[obss]
    y_test_1_vs_2 = df_test_1_vs_2['y']

    # Check for NaNs in the original data
    logging.info("Checking for NaNs in the original data...")
    logging.info(f"NaNs in X_train_0_vs_1: {np.isnan(X_train_0_vs_1).sum()}")
    logging.info(f"NaNs in X_test_0_vs_1: {np.isnan(X_test_0_vs_1).sum()}")
    logging.info(f"NaNs in X_train_0_vs_2: {np.isnan(X_train_0_vs_2).sum()}")
    logging.info(f"NaNs in X_test_0_vs_2: {np.isnan(X_test_0_vs_2).sum()}")
    logging.info(f"NaNs in X_train_1_vs_2: {np.isnan(X_train_1_vs_2).sum()}")
    logging.info(f"NaNs in X_test_1_vs_2: {np.isnan(X_test_1_vs_2).sum()}")
    logging.info(f"NaNs in X_test: {np.isnan(X_test).sum()}")

    # Impute missing values
    imputer = SimpleImputer(strategy='mean')
    X_train_0_vs_1 = imputer.fit_transform(X_train_0_vs_1)
    X_test_0_vs_1 = imputer.transform(X_test_0_vs_1)
    X_train_0_vs_2 = imputer.fit_transform(X_train_0_vs_2)
    X_test_0_vs_2 = imputer.transform(X_test_0_vs_2)
    X_train_1_vs_2 = imputer.fit_transform(X_train_1_vs_2)
    X_test_1_vs_2 = imputer.transform(X_test_1_vs_2)
    X_test = imputer.transform(X_test)

    # Check for remaining NaNs after imputation
    logging.info("Checking for NaNs after imputation...")
    logging.info(f"NaNs in X_train_0_vs_1: {np.isnan(X_train_0_vs_1).sum()}")
    logging.info(f"NaNs in X_test_0_vs_1: {np.isnan(X_test_0_vs_1).sum()}")
    logging.info(f"NaNs in X_train_0_vs_2: {np.isnan(X_train_0_vs_2).sum()}")
    logging.info(f"NaNs in X_test_0_vs_2: {np.isnan(X_test_0_vs_2).sum()}")
    logging.info(f"NaNs in X_train_1_vs_2: {np.isnan(X_train_1_vs_2).sum()}")
    logging.info(f"NaNs in X_test_1_vs_2: {np.isnan(X_test_1_vs_2).sum()}")
    logging.info(f"NaNs in X_test: {np.isnan(X_test).sum()}")

    # Set up the BDT models
    bdt_0_vs_1 = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
    bdt_0_vs_2 = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
    bdt_1_vs_2 = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)

    # Train the models
    logging.info("Training BDT models...")
    bdt_0_vs_1.fit(X_train_0_vs_1, y_train_0_vs_1)
    bdt_0_vs_2.fit(X_train_0_vs_2, y_train_0_vs_2)
    bdt_1_vs_2.fit(X_train_1_vs_2, y_train_1_vs_2)

    # Predict probabilities for the testing dataset
    y_pred_proba_0_vs_1 = bdt_0_vs_1.predict_proba(X_test)[:, 1]
    y_pred_proba_0_vs_2 = bdt_0_vs_2.predict_proba(X_test)[:, 1]
    y_pred_proba_1_vs_2 = bdt_1_vs_2.predict_proba(X_test)[:, 1]

    # Create a new DataFrame to store the actual y values and the three output scores
    df_scores = pd.DataFrame({
        'y_actual': y_test,
        'score_0_vs_1': y_pred_proba_0_vs_1,
        'score_0_vs_2': y_pred_proba_0_vs_2,
        'score_1_vs_2': y_pred_proba_1_vs_2
    })

    logging.info("Scores DataFrame head:")
    logging.info(df_scores.head())

    # Save the DataFrame to a CSV file
    df_scores.to_csv(f'{outputDir}/scores.csv', index=False)
    logging.info(f"Scores saved to {outputDir}/scores.csv")

    # Optional: Save the models for later use (if needed)
    joblib.dump(bdt_0_vs_1, f'{outputDir}/bdt_model_0_vs_1.pkl')
    joblib.dump(bdt_0_vs_2, f'{outputDir}/bdt_model_0_vs_2.pkl')
    joblib.dump(bdt_1_vs_2, f'{outputDir}/bdt_model_1_vs_2.pkl')
    logging.info("Models saved")

if __name__ == '__main__':
    main()