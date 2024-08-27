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
print("Parsing arguments")
args = parser.parse_args()
timeStamp = args.timestamp
outputDir = f'outputs/{timeStamp}'
coffeaFile = f"LHSskimmerOutput_{timeStamp}.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))


bbarIdx = 0
cbadIdx = 1
sbarIdx = 2
ubarIdx = 3
dbarIdx = 4
gIdx = 5
dIdx = 6
uIdx = 7
sIdx = 8
cIdx = 9
bIdx = 10

m_tt_edges = [400, 500, 600, 700, 900, 1200]
qqbarPairs = [(uIdx, ubarIdx), (dIdx, dbarIdx), (sIdx, sbarIdx), (cIdx, cbadIdx), (bIdx, bbarIdx)]


for era in out:
    eraHist = out[era]['yMatrix']

    # Initialize the A_FB matrix
    fbHist = eraHist.project('c', 'm_tt','p1', 'p2')

    # Loop through each bin in the 2D mesh
    A_FB = np.zeros(len(m_tt_edges) - 1)
    for i in range(len(m_tt_edges)-1):
        lowerL, higherL = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)
        cPos_qqbar = 0
        cPos_qbarq = 0
        cNeg_qqbar = 0
        cNeg_qbarq = 0
        for (qIdx, qbarIdx) in qqbarPairs:
            # print(i)
            cPos_qqbar += fbHist[0j:,lowerL:higherL,qIdx, qbarIdx].sum()
            cPos_qbarq += fbHist[0j:,lowerL:higherL,qbarIdx, qIdx].sum()
            cNeg_qqbar += fbHist[:0j,lowerL:higherL,qIdx, qbarIdx].sum() 
            cNeg_qbarq += fbHist[:0j,lowerL:higherL,qbarIdx, qIdx].sum()
        A_FB[i] = (((cPos_qqbar + cNeg_qbarq) - (cPos_qbarq + cNeg_qqbar))/((cPos_qqbar + cNeg_qbarq) + (cPos_qbarq + cNeg_qqbar)))

    print("The obtained A_FB values are: ", A_FB)

    cHist = eraHist.project('m_tt','deltay')

    A_C = np.zeros(len(m_tt_edges) - 1)
    for i in range(len(m_tt_edges)-1):
        lowerL, higherL = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)
        deltaPos = cHist[lowerL:higherL, 0j:].sum()
        deltaNeg = cHist[lowerL:higherL, :0j].sum()
        if (deltaPos + deltaNeg) != 0:
            A_C[i] = (deltaPos - deltaNeg)/(deltaPos + deltaNeg)
        else:
            A_C[i] = 0

    print("The obtained A_C values are: ", A_C)


    ytHist = eraHist.project('m_tt','y_t', 'y_tbar')


    A_out = np.zeros(len(m_tt_edges) - 1)
    A_in = np.zeros(len(m_tt_edges) - 1)
    y0 = 1.2j
    for i in range(len(m_tt_edges)-1):
        lowerL, higherL = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)
        yt_greaterThan_y0 = ytHist[lowerL:higherL, y0:, sum].sum()
        ytbar_greaterThan_y0 = ytHist[lowerL:higherL, sum, y0:].sum()
        yt_lowerThan_y0 = ytHist[lowerL:higherL, :y0, sum].sum()
        ytbar_lowerThan_y0 = ytHist[lowerL:higherL, sum, :y0].sum()

        if (yt_greaterThan_y0 + ytbar_greaterThan_y0) != 0:
            A_out[i] = (yt_greaterThan_y0 - ytbar_greaterThan_y0)/(yt_greaterThan_y0 + ytbar_greaterThan_y0)
        else:
            A_out[i] = 0


        if (yt_lowerThan_y0 + ytbar_lowerThan_y0) != 0:
            A_in[i] = (yt_lowerThan_y0 - ytbar_lowerThan_y0)/(yt_lowerThan_y0 + ytbar_lowerThan_y0)
        else:
            A_in[i] = 0

    print("The obtained A_out values are: ", A_out)

    print("The obtained A_in values are: ", A_in)


    np.save(f'{outputDir}/A_FB_{era}_{timeStamp}.npy', A_FB)
    np.save(f'{outputDir}/A_C_{era}_{timeStamp}.npy', A_C)
    np.save(f'{outputDir}/A_in_{era}_{timeStamp}.npy', A_in)
    np.save(f'{outputDir}/A_out_{era}_{timeStamp}.npy', A_out)


    np.save(f'{outputDir}/m_tt_edges_{era}_{timeStamp}.npy', m_tt_edges)


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