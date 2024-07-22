from coffea.util import load, save
import os, argparse
import hist
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description='Time stamp for book keeping')
parser.add_argument('-i', '--input', type=str, help='Specify the input filename', default='skimmerOutput.coffea')
parser.add_argument('-t', '--timestamp', type=str, help='Provide the time stamp for book keeping', default='tStamp')

args = parser.parse_args()

tStamp = args.timestamp


outputDir = 'outputs'
coffeaFile = args.input

# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

for era in ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']:
    print(f"Working on {era}")
    eraHist = out[era]['yMatrix']
    cHist = eraHist[:,:,:,:,:,:,0] ## Only qqbar channel
    proj = cHist.project('m_tt', 'beta_ttz')
    # Get the bin edges for m_tt and beta_ttz
    m_tt_edges = proj.axes[0].edges
    beta_ttz_edges = proj.axes[1].edges

    # Initialize the A_FB matrix
    A_FB = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    A_C = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    A_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    A_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    # Loop through each bin in the 2D mesh
    for i in range(len(m_tt_edges) - 1):
        for j in range(len(beta_ttz_edges) - 1):
            # Select the bin range for m_tt and beta_ttz
            m_tt_range = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)
            beta_ttz_range = (beta_ttz_edges[j]*1j, beta_ttz_edges[j + 1]*1j)
            counts_positive_c = cHist[0j:, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, sum].sum()
            counts_negative_c = cHist[:0j, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, sum].sum()
            counts_positive_delta = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, 0j:].sum()
            counts_negative_delta = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, :0j].sum()
            y0 = 1.2j
            counts_yt_greaterThan_y0 = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, sum, sum].sum()
            counts_yt_lowerThan_y0 = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, sum, sum].sum()
            counts_ytbar_greaterThan_y0 = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1],sum, y0:, sum].sum()
            counts_ytbar_lowerThan_y0 = cHist[sum, m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1],sum, :y0, sum].sum()


            # Calculate all the A's for this bin
            if (counts_positive_c + counts_negative_c) > 0:
                A_FB[i, j] = (counts_positive_c - counts_negative_c) / (counts_positive_c + counts_negative_c)
            else:
                A_FB[i, j] = 0.0  # Handle bins with zero total counts
                
            if (counts_positive_delta + counts_negative_delta) > 0:
                A_C[i, j] = (counts_positive_delta - counts_negative_delta) / (counts_positive_delta + counts_negative_delta)
            else:
                A_C[i, j] = 0.0  # Handle bins with zero total counts

            if (counts_yt_greaterThan_y0 + counts_ytbar_greaterThan_y0) > 0:
                A_out[i, j] = (counts_yt_greaterThan_y0 - counts_ytbar_greaterThan_y0) / (counts_yt_greaterThan_y0 + counts_ytbar_greaterThan_y0)
            else:
                A_out[i, j] = 0.0  # Handle bins with zero total counts


            if (counts_yt_lowerThan_y0 + counts_ytbar_lowerThan_y0) > 0:
                A_in[i, j] = (counts_yt_lowerThan_y0 - counts_ytbar_lowerThan_y0) / (counts_yt_lowerThan_y0 + counts_ytbar_lowerThan_y0)
            else:
                A_in[i, j] = 0.0  # Handle bins with zero total counts

    data = {
        'A_FB': '$A_{FB}$',
        'A_C': '$A_{C}$',
        'A_in': '$A_{in}$',
        'A_out': '$A_{out}$',
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
        plt.savefig(f'plots/{key}_{tStamp}.png')
        plt.close()

    AFB_flat = A_FB.flatten()
    AC_flat = A_C.flatten()
    Ain_flat = A_in.flatten()
    Aout_flat = A_out.flatten()

    np.save(f'outputs/AFBstar_flattened_array_{tStamp}.npy', AFB_flat)
    np.save(f'outputs/AC_flattened_array_{tStamp}.npy', AC_flat)
    np.save(f'outputs/Ain_flattened_array_{tStamp}.npy', Ain_flat)
    np.save(f'outputs/Aout_flattened_array_{tStamp}.npy', Aout_flat)

    dict = {
        'AC_flat': ('AC','$A_{{C}}$'),
        'Ain_flat': ('Ain','$A_{{in}}$'),
        'Aout_flat': ('Aout','$A_{{out}}$'),
    }

    # Loop through the dictionary and plot each array
    for key, (title, label) in dict.items():
        array = globals()[key]  # Access the array by name
        correlation, p_value = pearsonr(AFB_flat, array)
        print(f"\tFor {title}")
        print(f'\t\tCorrelation coefficient: {correlation}')
        print(f'\t\tP-value: {p_value}')
        # Optional: Plotting the data to visualize the correlation
        plt.figure()
        plt.scatter(AFB_flat, array)
        plt.xlabel('$A_{FB}$')
        plt.ylabel(label)
        plt.title(f'Correlation between $A_{{FB}}$ and {label}: {correlation:.2f}')
        plt.grid(True)
        plt.savefig(f'plots/AFB_{title}_{era}_{tStamp}.png')
        plt.close()