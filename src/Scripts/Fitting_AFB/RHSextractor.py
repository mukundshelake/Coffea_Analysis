from coffea.util import load, save
import os, argparse
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
coffeaFile = f"RHSskimmerOutput_{timeStamp}.coffea"


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

for era in out:
    fHist = out[era]['yMatrix']

    # Get the bin edges for m_tt and beta_ttz
    m_tt_edges = [300, 400, 500, 600, 700, 900, 1200]

    # Initialize the A_FB matrix
    Du = np.zeros(len(m_tt_edges) -1)
    Dd = np.zeros(len(m_tt_edges) -1)
    Fu_Nq = np.zeros(len(m_tt_edges) -1)
    Fd_Nq = np.zeros(len(m_tt_edges) -1)
    Fu_NEvents = np.zeros(len(m_tt_edges) -1)
    Fd_NEvents = np.zeros(len(m_tt_edges) -1)

    Du_out = np.zeros(len(m_tt_edges) -1)
    Dd_out = np.zeros(len(m_tt_edges) -1)
    Fu_Nq_out = np.zeros(len(m_tt_edges) -1)
    Fd_Nq_out = np.zeros(len(m_tt_edges) -1)
    Fu_NEvents_out = np.zeros(len(m_tt_edges) -1)
    Fd_NEvents_out = np.zeros(len(m_tt_edges) -1)
    outFactor = np.zeros(len(m_tt_edges) -1)


    Du_in = np.zeros(len(m_tt_edges) -1)
    Dd_in = np.zeros(len(m_tt_edges) -1)
    Fu_Nq_in = np.zeros(len(m_tt_edges) -1)
    Fd_Nq_in = np.zeros(len(m_tt_edges) -1)
    Fu_NEvents_in = np.zeros(len(m_tt_edges) -1)
    Fd_NEvents_in = np.zeros(len(m_tt_edges) -1)
    inFactor = np.zeros(len(m_tt_edges) -1)

    Du_notIN = np.zeros(len(m_tt_edges) -1)
    Dd_notIN = np.zeros(len(m_tt_edges) -1)
    Fu_Nq_notIN = np.zeros(len(m_tt_edges) -1)
    Fd_Nq_notIN = np.zeros(len(m_tt_edges) -1)
    Fu_NEvents_notIN = np.zeros(len(m_tt_edges) -1)
    Fd_NEvents_notIN = np.zeros(len(m_tt_edges) -1)


    Du_notOUT = np.zeros(len(m_tt_edges) -1)
    Dd_notOUT = np.zeros(len(m_tt_edges) -1)
    Fu_Nq_notOUT = np.zeros(len(m_tt_edges) -1)
    Fd_Nq_notOUT = np.zeros(len(m_tt_edges) -1)
    Fu_NEvents_notOUT = np.zeros(len(m_tt_edges) -1)
    Fd_NEvents_notOUT = np.zeros(len(m_tt_edges) -1)


    # Loop through each bin in the 2D mesh
    for i in range(len(m_tt_edges) - 1):
        # Select the bin range for m_tt and beta_ttz
        m_tt_range = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)

        Nuubar_highq1 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, uIdx, ubarIdx, 1].sum()
        Nuubar_highq2 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, uIdx, ubarIdx, 0].sum()

        Nubaru_highq1 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, ubarIdx, uIdx, 1].sum()
        Nubaru_highq2 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, ubarIdx, uIdx, 0].sum()

        Nuubar_highq = Nuubar_highq1 + Nubaru_highq2
        Nuubar_highqbar = Nuubar_highq2 + Nubaru_highq1


        Nddbar_highq1 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, dIdx, dbarIdx, 1].sum()
        Nddbar_highq2 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, dIdx, dbarIdx, 0].sum()

        Ndbard_highq1 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, dbarIdx, dIdx, 1].sum()
        Ndbard_highq2 = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, dbarIdx, dIdx, 0].sum()

        Nddbar_highq = Nddbar_highq1 + Ndbard_highq2
        Nddbar_highqbar = Nddbar_highq2 + Ndbard_highq1


        NEvents = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, sum, sum, sum].sum()
        N_gg = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, gIdx, gIdx, sum].sum()
        N_gx = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, gIdx, :, sum].sum()
        N_xg = fHist[m_tt_range[0]:m_tt_range[1], sum, sum, sum, :, gIdx, sum].sum()


        Nuubar = Nuubar_highq + Nuubar_highqbar
        Nddbar = Nddbar_highq + Nddbar_highqbar

        Nq = NEvents - N_gx - N_xg + N_gg


        if (Nuubar_highq + Nuubar_highqbar) > 0:
            Du[i] = (Nuubar_highq - Nuubar_highqbar)/(Nuubar_highq + Nuubar_highqbar)
        else:
            Du[i] = 0.0  # Handle bins with zero total counts

        if (Nddbar_highq + Nddbar_highqbar) > 0:
            Dd[i] = (Nddbar_highq - Nddbar_highqbar)/(Nddbar_highq + Nddbar_highqbar)
        else:
            Dd[i] = 0.0  # Handle bins with zero total counts

        if Nq > 0:
            Fu_Nq[i] = Nuubar/Nq
        else:
            Fu_Nq[i] = 0.0  # Handle bins with zero total counts

        if Nq > 0:
            Fd_Nq[i] = Nddbar/Nq
        else:
            Fd_Nq[i] = 0.0  # Handle bins with zero total counts

        if NEvents > 0:
            Fu_NEvents[i] = Nuubar/NEvents
        else:
            Fu_NEvents[i] = 0.0  # Handle bins with zero total counts

        if NEvents > 0:
            Fd_NEvents[i] = Nddbar/NEvents
        else:
            Fd_NEvents[i] = 0.0  # Handle bins with zero total counts


        y0 = 1.2j

        #### For region out
        Nuubar_highq1_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, uIdx,ubarIdx, 1].sum()
        Nuubar_highq2_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, uIdx,ubarIdx, 0].sum()

        Nubaru_highq1_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, ubarIdx,uIdx, 1].sum()
        Nubaru_highq2_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, ubarIdx,uIdx, 0].sum()

        Nuubar_highq_out = Nuubar_highq1_out + Nubaru_highq2_out
        Nuubar_highqbar_out = Nuubar_highq2_out + Nubaru_highq1_out


        Nddbar_highq1_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, dIdx,dbarIdx, 1].sum()
        Nddbar_highq2_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, dIdx,dbarIdx, 0].sum()

        Ndbard_highq1_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, dbarIdx,dIdx, 1].sum()
        Ndbard_highq2_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:,sum, dbarIdx,dIdx, 0].sum()

        Nddbar_highq_out = Nddbar_highq1_out + Ndbard_highq2_out
        Nddbar_highqbar_out = Nddbar_highq2_out + Ndbard_highq1_out



        NEvents_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, sum, sum, sum, sum].sum()
        N_gg_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, sum, gIdx, gIdx, sum].sum()
        N_gx_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, sum, gIdx, :, sum].sum()
        N_xg_out = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, sum, :, gIdx, sum].sum()


        Nuubar_out = Nuubar_highq_out + Nuubar_highqbar_out
        Nddbar_out = Nddbar_highq_out + Nddbar_highqbar_out

        Nq_out = NEvents_out - N_gx_out - N_xg_out + N_gg_out

        N_C = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, 0j:, sum, sum, sum].sum()
        N_Cprime = fHist[m_tt_range[0]:m_tt_range[1], y0:, y0:, :0j, sum, sum, sum].sum()
        N_B = fHist[m_tt_range[0]:m_tt_range[1], y0:, :y0, 0j:, sum, sum, sum].sum()
        N_D = fHist[m_tt_range[0]:m_tt_range[1], :y0, y0:, 0j:, sum, sum, sum].sum()


        if (Nuubar_highq_out + Nuubar_highqbar_out) > 0:
            Du_out[i] = (Nuubar_highq_out - Nuubar_highqbar_out)/(Nuubar_highq_out + Nuubar_highqbar_out)
        else:
            Du_out[i] = 0.0  # Handle bouts with zero total counts

        if (Nddbar_highq_out + Nddbar_highqbar_out) > 0:
            Dd_out[i] = (Nddbar_highq_out - Nddbar_highqbar_out)/(Nddbar_highq_out + Nddbar_highqbar_out)
        else:
            Dd_out[i] = 0.0  # Handle bouts with zero total counts

        if Nq_out > 0:
            Fu_Nq_out[i] = Nuubar_out/Nq_out
        else:
            Fu_Nq_out[i] = 0.0  # Handle bouts with zero total counts

        if Nq_out > 0:
            Fd_Nq_out[i] = Nddbar_out/Nq_out
        else:
            Fd_Nq_out[i] = 0.0  # Handle bouts with zero total counts

        if NEvents_out > 0:
            Fu_NEvents_out[i] = Nuubar_out/NEvents_out
        else:
            Fu_NEvents_out[i] = 0.0  # Handle bouts with zero total counts

        if NEvents_out > 0:
            Fd_NEvents_out[i] = Nddbar_out/NEvents_out
        else:
            Fd_NEvents_out[i] = 0.0  # Handle bouts with zero total counts

        if (N_B + N_D) > 0:
            outFactor[i] = 1/(1+2*((N_C + N_Cprime)/(N_B + N_D)))
        else:
            outFactor[i] = 0.0



        #### For region IN
        Nuubar_highq1_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, uIdx,ubarIdx, 1].sum()
        Nuubar_highq2_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, uIdx,ubarIdx, 0].sum()

        Nubaru_highq1_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, ubarIdx,uIdx, 1].sum()
        Nubaru_highq2_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, ubarIdx,uIdx, 0].sum()

        Nuubar_highq_in = Nuubar_highq1_in + Nubaru_highq2_in
        Nuubar_highqbar_in = Nuubar_highq2_in + Nubaru_highq1_in


        Nddbar_highq1_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, dIdx,dbarIdx, 1].sum()
        Nddbar_highq2_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, dIdx,dbarIdx, 0].sum()

        Ndbard_highq1_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, dbarIdx,dIdx, 1].sum()
        Ndbard_highq2_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0,sum, dbarIdx,dIdx, 0].sum()

        Nddbar_highq_in = Nddbar_highq1_in + Ndbard_highq2_in
        Nddbar_highqbar_in = Nddbar_highq2_in + Ndbard_highq1_in



        NEvents_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, sum, sum, sum, sum].sum()
        N_gg_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, sum, gIdx, gIdx, sum].sum()
        N_gx_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, sum, gIdx, :, sum].sum()
        N_xg_in = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, sum, :, gIdx, sum].sum()


        Nuubar_in = Nuubar_highq_in + Nuubar_highqbar_in
        Nddbar_in = Nddbar_highq_in + Nddbar_highqbar_in

        Nq_in = NEvents_in - N_gx_in - N_xg_in + N_gg_in

        N_A = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, 0j:, sum, sum, sum].sum()
        N_Aprime = fHist[m_tt_range[0]:m_tt_range[1], :y0, :y0, :0j, sum, sum, sum].sum()



        if (Nuubar_highq_in + Nuubar_highqbar_in) > 0:
            Du_in[i] = (Nuubar_highq_in - Nuubar_highqbar_in)/(Nuubar_highq_in + Nuubar_highqbar_in)
        else:
            Du_in[i] = 0.0  # Handle bouts with zero total counts

        if (Nddbar_highq_in + Nddbar_highqbar_in) > 0:
            Dd_in[i] = (Nddbar_highq_in - Nddbar_highqbar_in)/(Nddbar_highq_in + Nddbar_highqbar_in)
        else:
            Dd_in[i] = 0.0  # Handle bouts with zero total counts

        if Nq_in > 0:
            Fu_Nq_in[i] = Nuubar_in/Nq_in
        else:
            Fu_Nq_in[i] = 0.0  # Handle bouts with zero total counts

        if Nq_in > 0:
            Fd_Nq_in[i] = Nddbar_in/Nq_in
        else:
            Fd_Nq_in[i] = 0.0  # Handle bouts with zero total counts

        if NEvents_in > 0:
            Fu_NEvents_in[i] = Nuubar_in/NEvents_in
        else:
            Fu_NEvents_in[i] = 0.0  # Handle bouts with zero total counts

        if NEvents_in > 0:
            Fd_NEvents_in[i] = Nddbar_in/NEvents_in
        else:
            Fd_NEvents_in[i] = 0.0  # Handle bouts with zero total counts

        if (N_B + N_D) > 0:
            inFactor[i] = -1/(1+2*((N_A + N_Aprime)/(N_B + N_D)))
        else:
            inFactor[i] = 0.0
        
        #### For region notIN

        Nuubar_highq_notIN = Nuubar_highq - Nuubar_highq_in
        Nuubar_highqbar_notIN = Nuubar_highqbar - Nuubar_highqbar_in
        Nddbar_highq_notIN = Nddbar_highq - Nddbar_highq_in
        Nddbar_highqbar_notIN = Nddbar_highqbar - Nddbar_highqbar_in

        NEvents_notIN = NEvents - NEvents_in

        Nuubar_notIN = Nuubar_highq_notIN + Nuubar_highqbar_notIN
        Nddbar_notIN = Nddbar_highq_notIN + Nddbar_highqbar_notIN

        Nq_notIN = Nq - Nq_in

        if (Nuubar_highq_notIN + Nuubar_highqbar_notIN) > 0:
            Du_notIN[i] = (Nuubar_highq_notIN - Nuubar_highqbar_notIN)/(Nuubar_highq_notIN + Nuubar_highqbar_notIN)
        else:
            Du_notIN[i] = 0.0  # Handle bnotINs with zero total counts

        if (Nddbar_highq_notIN + Nddbar_highqbar_notIN) > 0:
            Dd_notIN[i] = (Nddbar_highq_notIN - Nddbar_highqbar_notIN)/(Nddbar_highq_notIN + Nddbar_highqbar_notIN)
        else:
            Dd_notIN[i] = 0.0  # Handle bnotINs with zero total counts

        if Nq_notIN > 0:
            Fu_Nq_notIN[i] = Nuubar_notIN/Nq_notIN
        else:
            Fu_Nq_notIN[i] = 0.0  # Handle bnotINs with zero total counts

        if Nq_notIN > 0:
            Fd_Nq_notIN[i] = Nddbar_notIN/Nq_notIN
        else:
            Fd_Nq_notIN[i] = 0.0  # Handle bnotINs with zero total counts

        if NEvents_notIN > 0:
            Fu_NEvents_notIN[i] = Nuubar_notIN/NEvents_notIN
        else:
            Fu_NEvents_notIN[i] = 0.0  # Handle bnotINs with zero total counts

        if NEvents_notIN > 0:
            Fd_NEvents_notIN[i] = Nddbar_notIN/NEvents_notIN
        else:
            Fd_NEvents_notIN[i] = 0.0  # Handle bnotINs with zero total counts


        #### For region notOUT
        Nuubar_highq_notOUT = Nuubar_highq - Nuubar_highq_out
        Nuubar_highqbar_notOUT = Nuubar_highqbar - Nuubar_highqbar_out
        Nddbar_highq_notOUT = Nddbar_highq - Nddbar_highq_out
        Nddbar_highqbar_notOUT = Nddbar_highqbar - Nddbar_highqbar_out

        NEvents_notOUT = NEvents - NEvents_out

        Nuubar_notOUT = Nuubar_highq_notOUT + Nuubar_highqbar_notOUT
        Nddbar_notOUT = Nddbar_highq_notOUT + Nddbar_highqbar_notOUT

        Nq_notOUT = Nq - Nq_out

        if (Nuubar_highq_notOUT + Nuubar_highqbar_notOUT) > 0:
            Du_notOUT[i] = (Nuubar_highq_notOUT - Nuubar_highqbar_notOUT)/(Nuubar_highq_notOUT + Nuubar_highqbar_notOUT)
        else:
            Du_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

        if (Nddbar_highq_notOUT + Nddbar_highqbar_notOUT) > 0:
            Dd_notOUT[i] = (Nddbar_highq_notOUT - Nddbar_highqbar_notOUT)/(Nddbar_highq_notOUT + Nddbar_highqbar_notOUT)
        else:
            Dd_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

        if Nq_notOUT > 0:
            Fu_Nq_notOUT[i] = Nuubar_notOUT/Nq_notOUT
        else:
            Fu_Nq_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

        if Nq_notOUT > 0:
            Fd_Nq_notOUT[i] = Nddbar_notOUT/Nq_notOUT
        else:
            Fd_Nq_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

        if NEvents_notOUT > 0:
            Fu_NEvents_notOUT[i] = Nuubar_notOUT/NEvents_notOUT
        else:
            Fu_NEvents_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

        if NEvents_notOUT > 0:
            Fd_NEvents_notOUT[i] = Nddbar_notOUT/NEvents_notOUT
        else:
            Fd_NEvents_notOUT[i] = 0.0  # Handle bnotOUTs with zero total counts

    np.save(f'{outputDir}/Du_{era}_{timeStamp}.npy', Du)
    np.save(f'{outputDir}/Dd_{era}_{timeStamp}.npy', Dd)
    np.save(f'{outputDir}/Fu_Nq_{era}_{timeStamp}.npy', Fu_Nq)
    np.save(f'{outputDir}/Fd_Nq_{era}_{timeStamp}.npy', Fd_Nq)
    np.save(f'{outputDir}/Fu_NEvents_{era}_{timeStamp}.npy', Fu_NEvents)
    np.save(f'{outputDir}/Fd_NEvents_{era}_{timeStamp}.npy', Fd_NEvents)

    np.save(f'{outputDir}/Du_out_{era}_{timeStamp}.npy', Du_out)
    np.save(f'{outputDir}/Dd_out_{era}_{timeStamp}.npy', Dd_out)
    np.save(f'{outputDir}/Fu_Nq_out_{era}_{timeStamp}.npy', Fu_Nq_out)
    np.save(f'{outputDir}/Fd_Nq_out_{era}_{timeStamp}.npy', Fd_Nq_out)
    np.save(f'{outputDir}/Fu_NEvents_out_{era}_{timeStamp}.npy', Fu_NEvents_out)
    np.save(f'{outputDir}/Fd_NEvents_out_{era}_{timeStamp}.npy', Fd_NEvents_out)
    np.save(f'{outputDir}/outFactor_{era}_{timeStamp}.npy', outFactor)

    np.save(f'{outputDir}/Du_in_{era}_{timeStamp}.npy', Du_in)
    np.save(f'{outputDir}/Dd_in_{era}_{timeStamp}.npy', Dd_in)
    np.save(f'{outputDir}/Fu_Nq_in_{era}_{timeStamp}.npy', Fu_Nq_in)
    np.save(f'{outputDir}/Fd_Nq_in_{era}_{timeStamp}.npy', Fd_Nq_in)
    np.save(f'{outputDir}/Fu_NEvents_in_{era}_{timeStamp}.npy', Fu_NEvents_in)
    np.save(f'{outputDir}/Fd_NEvents_in_{era}_{timeStamp}.npy', Fd_NEvents_in)
    np.save(f'{outputDir}/inFactor_{era}_{timeStamp}.npy', inFactor)

    np.save(f'{outputDir}/Du_notIN_{era}_{timeStamp}.npy', Du_notIN)
    np.save(f'{outputDir}/Dd_notIN_{era}_{timeStamp}.npy', Dd_notIN)
    np.save(f'{outputDir}/Fu_Nq_notIN_{era}_{timeStamp}.npy', Fu_Nq_notIN)
    np.save(f'{outputDir}/Fd_Nq_notIN_{era}_{timeStamp}.npy', Fd_Nq_notIN)
    np.save(f'{outputDir}/Fu_NEvents_notIN_{era}_{timeStamp}.npy', Fu_NEvents_notIN)
    np.save(f'{outputDir}/Fd_NEvents_notIN_{era}_{timeStamp}.npy', Fd_NEvents_notIN)


    np.save(f'{outputDir}/Du_notOUT_{era}_{timeStamp}.npy', Du_notOUT)
    np.save(f'{outputDir}/Dd_notOUT_{era}_{timeStamp}.npy', Dd_notOUT)
    np.save(f'{outputDir}/Fu_Nq_notOUT_{era}_{timeStamp}.npy', Fu_Nq_notOUT)
    np.save(f'{outputDir}/Fd_Nq_notOUT_{era}_{timeStamp}.npy', Fd_Nq_notOUT)
    np.save(f'{outputDir}/Fu_NEvents_notOUT_{era}_{timeStamp}.npy', Fu_NEvents_notOUT)
    np.save(f'{outputDir}/Fd_NEvents_notOUT_{era}_{timeStamp}.npy', Fd_NEvents_notOUT)

    np.save(f'{outputDir}/m_tt_edges_{era}_{timeStamp}.npy', m_tt_edges)

    # Flatten the 2D arrays
    Du_flat = Du.flatten()
    Dd_flat = Dd.flatten()
    Fu_Nq_flat = Fu_Nq.flatten()
    Fd_Nq_flat = Fd_Nq.flatten()
    Fu_NEvents_flat = Fu_NEvents.flatten()
    Fd_NEvents_flat = Fd_NEvents.flatten()

    Du_out_flat = Du_out.flatten()
    Dd_out_flat = Dd_out.flatten()
    Fu_Nq_out_flat = Fu_Nq_out.flatten()
    Fd_Nq_out_flat = Fd_Nq_out.flatten()
    Fu_NEvents_out_flat = Fu_NEvents_out.flatten()
    Fd_NEvents_out_flat = Fd_NEvents_out.flatten()
    outFactor_flat = outFactor.flatten()

    Du_in_flat = Du_in.flatten()
    Dd_in_flat = Dd_in.flatten()
    Fu_Nq_in_flat = Fu_Nq_in.flatten()
    Fd_Nq_in_flat = Fd_Nq_in.flatten()
    Fu_NEvents_in_flat = Fu_NEvents_in.flatten()
    Fd_NEvents_in_flat = Fd_NEvents_in.flatten()
    inFactor_flat = inFactor.flatten()

    Du_notIN_flat = Du_notIN.flatten()
    Dd_notIN_flat = Dd_notIN.flatten()
    Fu_Nq_notIN_flat = Fu_Nq_notIN.flatten()
    Fd_Nq_notIN_flat = Fd_Nq_notIN.flatten()
    Fu_NEvents_notIN_flat = Fu_NEvents_notIN.flatten()
    Fd_NEvents_notIN_flat = Fd_NEvents_notIN.flatten()

    Du_notOUT_flat = Du_notOUT.flatten()
    Dd_notOUT_flat = Dd_notOUT.flatten()
    Fu_Nq_notOUT_flat = Fu_Nq_notOUT.flatten()
    Fd_Nq_notOUT_flat = Fd_Nq_notOUT.flatten()
    Fu_NEvents_notOUT_flat = Fu_NEvents_notOUT.flatten()
    Fd_NEvents_notOUT_flat = Fd_NEvents_notOUT.flatten()

    # Create a DataFrame
    df = pd.DataFrame({
        'Du': Du_flat,
        'Dd': Dd_flat,
        'Fu_Nq': Fu_Nq_flat,
        'Fd_Nq': Fd_Nq_flat,
        'Fu_NEvents': Fu_NEvents_flat,
        'Fd_NEvents': Fd_NEvents_flat,
        'Du_out': Du_out_flat,
        'Dd_out': Dd_out_flat,
        'Fu_Nq_out': Fu_Nq_out_flat,
        'Fd_Nq_out': Fd_Nq_out_flat,
        'Fu_NEvents_out': Fu_NEvents_out_flat,
        'Fd_NEvents_out': Fd_NEvents_out_flat,
        'outFactor': outFactor_flat,
        'Du_in': Du_in_flat,
        'Dd_in': Dd_in_flat,
        'Fu_Nq_in': Fu_Nq_in_flat,
        'Fd_Nq_in': Fd_Nq_in_flat,
        'Fu_NEvents_in': Fu_NEvents_in_flat,
        'Fd_NEvents_in': Fd_NEvents_in_flat,
        'inFactor': inFactor_flat,
        'Du_notIN': Du_notIN_flat,
        'Dd_notIN': Dd_notIN_flat,
        'Fu_Nq_notIN': Fu_Nq_notIN_flat,
        'Fd_Nq_notIN': Fd_Nq_notIN_flat,
        'Fu_NEvents_notIN': Fu_NEvents_notIN_flat,
        'Fd_NEvents_notIN': Fd_NEvents_notIN_flat,
        'Du_notOUT': Du_notOUT_flat,
        'Dd_notOUT': Dd_notOUT_flat,
        'Fu_Nq_notOUT': Fu_Nq_notOUT_flat,
        'Fd_Nq_notOUT': Fd_Nq_notOUT_flat,
        'Fu_NEvents_notOUT': Fu_NEvents_notOUT_flat,
        'Fd_NEvents_notOUT': Fd_NEvents_notOUT_flat
    })

    # Print DataFrame info and first few rows
    print(df.info())
    print(df.head())


    # Define the filename for the CSV file
    csv_filename = f'{outputDir}/RHSdata_{era}_{timeStamp}.csv'  # Customize as needed

    # Save the DataFrame to a CSV file
    df.to_csv(csv_filename, index=False)

    print(f"DataFrame saved to {csv_filename}")