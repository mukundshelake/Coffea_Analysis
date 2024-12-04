import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import logging
import os

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
    parser = argparse.ArgumentParser(description="Plot results from BDT scores.")
    parser.add_argument(
        '-t', '--timestamp',
        type=str,
        required=True,
        help="Specify the timestamp used in the scores CSV file."
    )
    args = parser.parse_args()

    timestamp = args.timestamp
    outputDir = f'outputs/{timestamp}'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    setup_logging(os.path.splitext(os.path.basename(__file__))[0], outputDir)

    logging.info(f"Reading scores from {outputDir}/scores.csv")
    df_scores = pd.read_csv(f'{outputDir}/scores.csv')

    # Rename the columns for clarity
    df_scores.rename(columns={
        'score_0_vs_1': 'N01',
        'score_0_vs_2': 'N02',
        'score_1_vs_2': 'N12'
    }, inplace=True)

    # Create new columns for the differences
    df_scores['N01_N02'] = df_scores['N01'] - df_scores['N02']
    df_scores['N02_N12'] = df_scores['N02'] - df_scores['N12']
    df_scores['N12_N01'] = df_scores['N12'] - df_scores['N01']

    # Define pairs to plot
    pairs = [
        ('N01', 'N02'),
        ('N01', 'N12'),
        ('N02', 'N12'),
        ('N01_N02', 'N02_N12'),
        ('N02_N12', 'N12_N01'),
        ('N12_N01', 'N01_N02')
    ]

    logging.info(f"Number of pairs to plot: {len(pairs)}")

    # Plot each pair
    for x, y in pairs:
        logging.info(f"Plotting pair: {x} vs {y}")
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # All y_actual values
        sns.scatterplot(data=df_scores, x=x, y=y, hue='y_actual', palette='viridis', alpha=0.6, ax=axes[0, 0])
        axes[0, 0].set_title(f'Scatter Plot: {x} vs {y} (All)')
        
        # y_actual == 0
        sns.scatterplot(data=df_scores[df_scores['y_actual'] == 0], x=x, y=y, palette='viridis', alpha=0.6, ax=axes[0, 1])
        axes[0, 1].set_title(f'Scatter Plot: {x} vs {y} (y_actual == 0)')
        
        # y_actual == 1
        sns.scatterplot(data=df_scores[df_scores['y_actual'] == 1], x=x, y=y, palette='viridis', alpha=0.6, ax=axes[1, 0])
        axes[1, 0].set_title(f'Scatter Plot: {x} vs {y} (y_actual == 1)')
        
        # y_actual == 2
        sns.scatterplot(data=df_scores[df_scores['y_actual'] == 2], x=x, y=y, palette='viridis', alpha=0.6, ax=axes[1, 1])
        axes[1, 1].set_title(f'Scatter Plot: {x} vs {y} (y_actual == 2)')
        
        plt.tight_layout()
        plt.savefig(f'{outputDir}/scatter_plot_grid_{x}_vs_{y}.png')
        logging.info(f"Scatter plot saved to {outputDir}/scatter_plot_grid_{x}_vs_{y}.png")
        plt.clf()

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # All y_actual values
        sns.kdeplot(data=df_scores, x=x, y=y, levels=5, color='k', linewidths=1, ax=axes[0, 0])
        axes[0, 0].set_title(f'Density Plot: {x} vs {y} (All)')
        
        # y_actual == 0
        sns.kdeplot(data=df_scores[df_scores['y_actual'] == 0], x=x, y=y, levels=5, color='k', linewidths=1, ax=axes[0, 1])
        axes[0, 1].set_title(f'Density Plot: {x} vs {y} (y_actual == 0)')
        
        # y_actual == 1
        sns.kdeplot(data=df_scores[df_scores['y_actual'] == 1], x=x, y=y, levels=5, color='k', linewidths=1, ax=axes[1, 0])
        axes[1, 0].set_title(f'Density Plot: {x} vs {y} (y_actual == 1)')
        
        # y_actual == 2
        sns.kdeplot(data=df_scores[df_scores['y_actual'] == 2], x=x, y=y, levels=5, color='k', linewidths=1, ax=axes[1, 1])
        axes[1, 1].set_title(f'Density Plot: {x} vs {y} (y_actual == 2)')
        
        plt.tight_layout()
        plt.savefig(f'{outputDir}/density_plot_grid_{x}_vs_{y}.png')
        logging.info(f"Density plot saved to {outputDir}/density_plot_grid_{x}_vs_{y}.png")
        plt.clf()

if __name__ == '__main__':
    main()