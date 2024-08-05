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
outputDir = 'outputs'
coffeaFile = f"skimmerOutput_{timeStamp}.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

for era in out:
    # Load the .npy files
    Du = np.load(f'{outputDir}/Du_{era}_{timeStamp}.npy')
    Dd = np.load(f'{outputDir}/Dd_{era}_{timeStamp}.npy')
    Fu_Nq = np.load(f'{outputDir}/Fu_Nq_{era}_{timeStamp}.npy')
    Fd_Nq = np.load(f'{outputDir}/Fd_Nq_{era}_{timeStamp}.npy')
    Fu_NEvents = np.load(f'{outputDir}/Fu_NEvents_{era}_{timeStamp}.npy')
    Fd_NEvents = np.load(f'{outputDir}/Fd_NEvents_{era}_{timeStamp}.npy')

    Du_out = np.load(f'{outputDir}/Du_out_{era}_{timeStamp}.npy')
    Dd_out = np.load(f'{outputDir}/Dd_out_{era}_{timeStamp}.npy')
    Fu_Nq_out = np.load(f'{outputDir}/Fu_Nq_out_{era}_{timeStamp}.npy')
    Fd_Nq_out = np.load(f'{outputDir}/Fd_Nq_out_{era}_{timeStamp}.npy')
    Fu_NEvents_out = np.load(f'{outputDir}/Fu_NEvents_out_{era}_{timeStamp}.npy')
    Fd_NEvents_out = np.load(f'{outputDir}/Fd_NEvents_out_{era}_{timeStamp}.npy')

    Du_in = np.load(f'{outputDir}/Du_in_{era}_{timeStamp}.npy')
    Dd_in = np.load(f'{outputDir}/Dd_in_{era}_{timeStamp}.npy')
    Fu_Nq_in = np.load(f'{outputDir}/Fu_Nq_in_{era}_{timeStamp}.npy')
    Fd_Nq_in = np.load(f'{outputDir}/Fd_Nq_in_{era}_{timeStamp}.npy')
    Fu_NEvents_in = np.load(f'{outputDir}/Fu_NEvents_in_{era}_{timeStamp}.npy')
    Fd_NEvents_in = np.load(f'{outputDir}/Fd_NEvents_in_{era}_{timeStamp}.npy')

    Du_notIN = np.load(f'{outputDir}/Du_notIN_{era}_{timeStamp}.npy')
    Dd_notIN = np.load(f'{outputDir}/Dd_notIN_{era}_{timeStamp}.npy')
    Fu_Nq_notIN = np.load(f'{outputDir}/Fu_Nq_notIN_{era}_{timeStamp}.npy')
    Fd_Nq_notIN = np.load(f'{outputDir}/Fd_Nq_notIN_{era}_{timeStamp}.npy')
    Fu_NEvents_notIN = np.load(f'{outputDir}/Fu_NEvents_notIN_{era}_{timeStamp}.npy')
    Fd_NEvents_notIN = np.load(f'{outputDir}/Fd_NEvents_notIN_{era}_{timeStamp}.npy')

    Du_notOUT = np.load(f'{outputDir}/Du_notOUT_{era}_{timeStamp}.npy')
    Dd_notOUT = np.load(f'{outputDir}/Dd_notOUT_{era}_{timeStamp}.npy')
    Fu_Nq_notOUT = np.load(f'{outputDir}/Fu_Nq_notOUT_{era}_{timeStamp}.npy')
    Fd_Nq_notOUT = np.load(f'{outputDir}/Fd_Nq_notOUT_{era}_{timeStamp}.npy')
    Fu_NEvents_notOUT = np.load(f'{outputDir}/Fu_NEvents_notOUT_{era}_{timeStamp}.npy')
    Fd_NEvents_notOUT = np.load(f'{outputDir}/Fd_NEvents_notOUT_{era}_{timeStamp}.npy')

    m_tt_edges = np.load(f'{outputDir}/m_tt_edges_{era}_{timeStamp}.npy')

    beta_ttz_edges = np.load(f'{outputDir}/beta_ttz_edges_{era}_{timeStamp}.npy')


    data = {
        'Du': '$D_{u}$',
        'Dd': '$D_{d}$',
        'Fu_Nq': '$F_{u}\ w.r.t.\ N_{q}$',
        'Fd_Nq': '$F_{d}\ w.r.t.\ N_{q}$',
        'Fu_NEvents': '$F_{u}\ w.r.t.\ N_{total}$',
        'Fd_NEvents': '$F_{d}\ w.r.t.\ N_{total}$',
        'Du_out': '$Out\ Region\ D_{u}$',
        'Dd_out': '$Out\ Region\ D_{d}$',
        'Fu_Nq_out': '$Out\ Region\ F_{u}\ w.r.t.\ N_{q}$',
        'Fd_Nq_out': '$Out\ Region\ F_{d}\ w.r.t.\ N_{q}$',
        'Fu_NEvents_out': '$Out\ Region\ F_{u}\ w.r.t.\ N_{total}$',
        'Fd_NEvents_out': '$Out\ Region\ F_{d}\ w.r.t.\ N_{total}$',
        'Du_in': '$In\ Region\ D_{u}$',
        'Dd_in': '$In\ Region\ D_{d}$',
        'Fu_Nq_in': '$In\ Region\ F_{u}\ w.r.t.\ N_{q}$',
        'Fd_Nq_in': '$In\ Region\ F_{d}\ w.r.t.\ N_{q}$',
        'Fu_NEvents_in': '$In\ Region\ F_{u}\ w.r.t.\ N_{total}$',
        'Fd_NEvents_in': '$In\ Region\ F_{d}\ w.r.t.\ N_{total}$',
        'Du_notIN': '$Non-In\ Region\ D_{u}$',
        'Dd_notIN': '$Non-In\ Region\ D_{d}$',
        'Fu_Nq_notIN': '$Non-In\ Region\ F_{u}\ w.r.t.\ N_{q}$',
        'Fd_Nq_notIN': '$Non-In\ Region\ F_{d}\ w.r.t.\ N_{q}$',
        'Fu_NEvents_notIN': '$Non-In\ Region\ F_{u}\ w.r.t.\ N_{total}$',
        'Fd_NEvents_notIN': '$Non-In\ Region\ F_{d}\ w.r.t.\ N_{total}$',
        'Du_notOUT': '$Non-Out\ Region\ D_{u}$',
        'Dd_notOUT': '$Non-Out\ Region\ D_{d}$',
        'Fu_Nq_notOUT': '$Non-Out\ Region\ F_{u}\ w.r.t.\ N_{q}$',
        'Fd_Nq_notOUT': '$Non-Out\ Region\ F_{d}\ w.r.t.\ N_{q}$',
        'Fu_NEvents_notOUT': '$Non-Out\ Region\ F_{u}\ w.r.t.\ N_{total}$',
        'Fd_NEvents_notOUT': '$Non-Out\ Region\ F_{d}\ w.r.t.\ N_{total}$',
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