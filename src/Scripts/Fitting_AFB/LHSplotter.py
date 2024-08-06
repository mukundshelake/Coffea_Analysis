from coffea.util import load, save
import os, argparse
import hist
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description="Process some eras.")

parser.add_argument(
    '-t','--timestamp',
    type=str,
    default='timestamp',
    help="Specify the timestamp. Default is 'timestamp'."
)

args = parser.parse_args()
timeStamp = args.timestamp
outputDir = f'outputs/{timeStamp}'
coffeaFile = f"LHSskimmerOutput_{timeStamp}.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

for era in out:
    A_FB = np.load(f'{outputDir}/A_FB_{era}_{timeStamp}.npy')
    A_C = np.load(f'{outputDir}/A_C_{era}_{timeStamp}.npy')
    A_in = np.load(f'{outputDir}/A_in_{era}_{timeStamp}.npy')
    A_out = np.load(f'{outputDir}/A_out_{era}_{timeStamp}.npy')

    m_tt_edges = np.load(f'{outputDir}/m_tt_edges_{era}_{timeStamp}.npy')

    beta_ttz_edges = np.load(f'{outputDir}/beta_ttz_edges_{era}_{timeStamp}.npy')


    data = {
    'A_FB': '$A_{FB}$',
    'A_C': '$A_{C}$',
    'A_in': '$A_{in}$',
    'A_out': '$A_{out}$'
    }

    # Loop through the dictionary and plot each array
    for key, value in data.items():
        array = globals()[key]  # Access the array by name
        X, Y = np.meshgrid(m_tt_edges[:-1], beta_ttz_edges[:-1])
        fig, ax = plt.subplots()
        c = ax.pcolormesh(X, Y, array.T, cmap='coolwarm', shading='auto')
        ax.set_xlabel('$m_{tt}$')
        ax.set_ylabel('$beta_{tt}^z$')
        plt.title(value)
        fig.colorbar(c, ax=ax, label=value)
        plt.savefig(f'results/plots/{key}_{era}_{timeStamp}.png')
        # plt.show()
        plt.close()