from coffea.util import load, save
import os, argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    cHist = out[era]['yMatrix']

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


    np.save(f'{outputDir}/A_FB_{era}_{timeStamp}.npy', A_FB)
    np.save(f'{outputDir}/A_C_{era}_{timeStamp}.npy', A_C)
    np.save(f'{outputDir}/A_in_{era}_{timeStamp}.npy', A_in)
    np.save(f'{outputDir}/A_out_{era}_{timeStamp}.npy', A_out)


    np.save(f'{outputDir}/m_tt_edges_{era}_{timeStamp}.npy', m_tt_edges)
    np.save(f'{outputDir}/beta_ttz_edges_{era}_{timeStamp}.npy', beta_ttz_edges)


    A_FB_flat = A_FB.flatten()
    A_C_flat = A_C.flatten()
    A_in_flat = A_in.flatten()
    A_out_flat = A_out.flatten()


    df = pd.DataFrame({
        'A_FB': A_FB_flat,
        'A_C': A_C_flat,
        'A_in': A_in_flat,
        'A_out': A_out_flat
    })
    # Print DataFrame info and first few rows
    print(df.info())
    print(df.head())

    # Define the filename for the CSV file
    csv_filename = f'{outputDir}/LHSdata_{era}_{timeStamp}.csv'  # Customize as needed

    # Save the DataFrame to a CSV file
    df.to_csv(csv_filename, index=False)

    print(f"DataFrame saved to {csv_filename}")