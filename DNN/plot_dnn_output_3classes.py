import json
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_dnn_output(json_file, output_dir):
    # Load the output data from the JSON file
    with open(json_file, 'r') as f:
        output_data = json.load(f)

    # Extract y_test and y_pred
    y_test = output_data['y_test']
    y_pred = output_data['y_pred']

    # Convert y_pred to a DataFrame
    y_pred_df = pd.DataFrame(y_pred, columns=['N1', 'N2', 'N3'])
    y_pred_df['y_test'] = y_test

    # Calculate differences
    y_pred_df['N1-N2'] = y_pred_df['N1'] - y_pred_df['N2']
    y_pred_df['N2-N3'] = y_pred_df['N2'] - y_pred_df['N3']
    y_pred_df['N3-N1'] = y_pred_df['N3'] - y_pred_df['N1']

    # Define plot pairs
    plot_pairs = [
        ('N1', 'N2'),
        ('N1', 'N3'),
        ('N2', 'N3'),
        ('N1-N2', 'N2-N3'),
        ('N2-N3', 'N3-N1'),
        ('N3-N1', 'N1-N2')
    ]

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Define class name mapping
    class_name_map = {0: 'qq', 1: 'gg', 2: 'qg'}

    # Create scatter plots
    for x, y in plot_pairs:
        fig, axs = plt.subplots(2, 2, figsize=(15, 12))

        # Integrated scatter plot
        sns.scatterplot(ax=axs[0, 0], x=x, y=y, hue='y_test', palette='viridis', data=y_pred_df, alpha=0.5)
        axs[0, 0].set_title(f'Integrated scatter plot of {x} vs {y} with y_test as color')
        axs[0, 0].set_xlabel(x)
        axs[0, 0].set_ylabel(y)
        axs[0, 0].legend(title='Actual Class')

        # Separate scatter plots for each class
        for class_value, ax in zip([0, 1, 2], [axs[0, 1], axs[1, 0], axs[1, 1]]):
            sns.scatterplot(ax=ax, x=x, y=y, hue='y_test', palette='viridis', data=y_pred_df[y_pred_df['y_test'] == class_value], alpha=0.5)
            ax.set_title(f'Scatter plot of {x} vs {y} for class {class_name_map[class_value]}')
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.legend(title='Actual Class')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'scatter_plot_{x}_vs_{y}.png'))
        plt.show()
        plt.close()

    # Create density plots
    for x, y in plot_pairs:
        fig, axs = plt.subplots(2, 2, figsize=(15, 12))

        # Integrated density plot
        sns.kdeplot(ax=axs[0, 0], x=x, y=y, hue='y_test', palette='viridis', data=y_pred_df, fill=True, alpha=0.3)
        axs[0, 0].set_title(f'Integrated density plot of {x} vs {y} with y_test as color')
        axs[0, 0].set_xlabel(x)
        axs[0, 0].set_ylabel(y)
        handles, labels = axs[0, 0].get_legend_handles_labels()
        axs[0, 0].legend(handles=handles, labels=['Class qq', 'Class gg', 'Class qg'], title='Actual Class')

        # Separate density plots for each class
        for class_value, ax in zip([0, 1, 2], [axs[0, 1], axs[1, 0], axs[1, 1]]):
            sns.kdeplot(ax=ax, x=x, y=y, hue='y_test', palette='viridis', data=y_pred_df[y_pred_df['y_test'] == class_value], fill=True, alpha=0.3)
            ax.set_title(f'Density plot of {x} vs {y} for class {class_name_map[class_value]}')
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles=handles, labels=[f'Class {class_name_map[class_value]}'], title='Actual Class')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'density_plot_{x}_vs_{y}.png'))
        plt.show()
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot DNN output.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the JSON file containing the DNN output.')
    parser.add_argument('-o', '--output', type=str, default='outputs', help='Output directory for the plots.')
    parser.add_argument('-t', '--timestamp', type=str, default='timestamp', help='Timestamp to create a subfolder for the plots.')

    args = parser.parse_args()

    output_dir = os.path.join(args.output, args.timestamp)
    plot_dnn_output(args.file, output_dir)
