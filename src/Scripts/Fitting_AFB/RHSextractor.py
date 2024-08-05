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
outputDir = 'outputs'
coffeaFile = f"skimmerOutput_{timeStamp}.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))


for era in out:
    eraHist = out[era]['yMatrix']
    fHist = eraHist[sum,:,:,:,:,sum,:]

    proj = fHist.project('m_tt', 'beta_ttz')
    # Get the bin edges for m_tt and beta_ttz
    m_tt_edges = proj.axes[0].edges
    beta_ttz_edges = proj.axes[1].edges

    # Initialize the A_FB matrix
    Du = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Dd = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_Nq = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_Nq = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_NEvents = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_NEvents = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))

    Du_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Dd_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_Nq_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_Nq_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_NEvents_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_NEvents_out = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))


    Du_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Dd_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_Nq_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_Nq_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_NEvents_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_NEvents_in = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))

    Du_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Dd_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_Nq_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_Nq_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_NEvents_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_NEvents_notIN = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))


    Du_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Dd_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_Nq_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_Nq_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fu_NEvents_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))
    Fd_NEvents_notOUT = np.zeros((len(m_tt_edges) - 1, len(beta_ttz_edges) - 1))


    # Loop through each bin in the 2D mesh
    for i in range(len(m_tt_edges) - 1):
        for j in range(len(beta_ttz_edges) - 1):
            # Select the bin range for m_tt and beta_ttz
            m_tt_range = (m_tt_edges[i]*1j, m_tt_edges[i + 1]*1j)
            beta_ttz_range = (beta_ttz_edges[j]*1j, beta_ttz_edges[j + 1]*1j)
            Nuubar_highq = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, 1].sum()
            Nuubar_highqbar = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, -1].sum()
            Nddbar_highq = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, 2].sum()
            Nddbar_highqbar = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, -2].sum()

            NEvents = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, sum].sum()
            N_non_udQ = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], sum, sum, 0].sum()

            Nuubar = Nuubar_highq + Nuubar_highqbar
            Nddbar = Nddbar_highq + Nddbar_highqbar

            Nq = Nuubar + Nddbar + N_non_udQ


            if (Nuubar_highq + Nuubar_highqbar) > 0:
                Du[i, j] = (Nuubar_highq - Nuubar_highqbar)/(Nuubar_highq + Nuubar_highqbar)
            else:
                Du[i, j] = 0.0  # Handle bins with zero total counts

            if (Nddbar_highq + Nddbar_highqbar) > 0:
                Dd[i, j] = (Nddbar_highq - Nddbar_highqbar)/(Nddbar_highq + Nddbar_highqbar)
            else:
                Dd[i, j] = 0.0  # Handle bins with zero total counts

            if Nq > 0:
                Fu_Nq[i, j] = Nuubar/Nq
            else:
                Fu_Nq[i, j] = 0.0  # Handle bins with zero total counts

            if Nq > 0:
                Fd_Nq[i, j] = Nddbar/Nq
            else:
                Fd_Nq[i, j] = 0.0  # Handle bins with zero total counts

            if NEvents > 0:
                Fu_NEvents[i, j] = Nuubar/NEvents
            else:
                Fu_NEvents[i, j] = 0.0  # Handle bins with zero total counts

            if NEvents > 0:
                Fd_NEvents[i, j] = Nddbar/NEvents
            else:
                Fd_NEvents[i, j] = 0.0  # Handle bins with zero total counts


            y0 = 1.2j

            #### For region out
            Nuubar_highq_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, 1].sum()
            Nuubar_highqbar_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, -1].sum()
            Nddbar_highq_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, 2].sum()
            Nddbar_highqbar_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, -2].sum()

            NEvents_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, sum].sum()
            N_non_udQ_out = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], y0:, y0:, 0].sum()

            Nuubar_out = Nuubar_highq_out + Nuubar_highqbar_out
            Nddbar_out = Nddbar_highq_out + Nddbar_highqbar_out

            Nq_out = Nuubar_out + Nddbar_out + N_non_udQ_out

            if (Nuubar_highq_out + Nuubar_highqbar_out) > 0:
                Du_out[i, j] = (Nuubar_highq_out - Nuubar_highqbar_out)/(Nuubar_highq_out + Nuubar_highqbar_out)
            else:
                Du_out[i, j] = 0.0  # Handle bouts with zero total counts

            if (Nddbar_highq_out + Nddbar_highqbar_out) > 0:
                Dd_out[i, j] = (Nddbar_highq_out - Nddbar_highqbar_out)/(Nddbar_highq_out + Nddbar_highqbar_out)
            else:
                Dd_out[i, j] = 0.0  # Handle bouts with zero total counts

            if Nq_out > 0:
                Fu_Nq_out[i, j] = Nuubar_out/Nq_out
            else:
                Fu_Nq_out[i, j] = 0.0  # Handle bouts with zero total counts

            if Nq_out > 0:
                Fd_Nq_out[i, j] = Nddbar_out/Nq_out
            else:
                Fd_Nq_out[i, j] = 0.0  # Handle bouts with zero total counts

            if NEvents_out > 0:
                Fu_NEvents_out[i, j] = Nuubar_out/NEvents_out
            else:
                Fu_NEvents_out[i, j] = 0.0  # Handle bouts with zero total counts

            if NEvents_out > 0:
                Fd_NEvents_out[i, j] = Nddbar_out/NEvents_out
            else:
                Fd_NEvents_out[i, j] = 0.0  # Handle bouts with zero total counts


            #### For region IN

            Nuubar_highq_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, 1].sum()
            Nuubar_highqbar_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, -1].sum()
            Nddbar_highq_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, 2].sum()
            Nddbar_highqbar_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, -2].sum()

            NEvents_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, sum].sum()
            N_non_udQ_in = fHist[m_tt_range[0]:m_tt_range[1], beta_ttz_range[0]: beta_ttz_range[1], :y0, :y0, 0].sum()

            Nuubar_in = Nuubar_highq_in + Nuubar_highqbar_in
            Nddbar_in = Nddbar_highq_in + Nddbar_highqbar_in

            Nq_in = Nuubar_in + Nddbar_in + N_non_udQ_in

            if (Nuubar_highq_in + Nuubar_highqbar_in) > 0:
                Du_in[i, j] = (Nuubar_highq_in - Nuubar_highqbar_in)/(Nuubar_highq_in + Nuubar_highqbar_in)
            else:
                Du_in[i, j] = 0.0  # Handle bins with zero total counts

            if (Nddbar_highq_in + Nddbar_highqbar_in) > 0:
                Dd_in[i, j] = (Nddbar_highq_in - Nddbar_highqbar_in)/(Nddbar_highq_in + Nddbar_highqbar_in)
            else:
                Dd_in[i, j] = 0.0  # Handle bins with zero total counts

            if Nq_in > 0:
                Fu_Nq_in[i, j] = Nuubar_in/Nq_in
            else:
                Fu_Nq_in[i, j] = 0.0  # Handle bins with zero total counts

            if Nq_in > 0:
                Fd_Nq_in[i, j] = Nddbar_in/Nq_in
            else:
                Fd_Nq_in[i, j] = 0.0  # Handle bins with zero total counts

            if NEvents_in > 0:
                Fu_NEvents_in[i, j] = Nuubar_in/NEvents_in
            else:
                Fu_NEvents_in[i, j] = 0.0  # Handle bins with zero total counts

            if NEvents_in > 0:
                Fd_NEvents_in[i, j] = Nddbar_in/NEvents_in
            else:
                Fd_NEvents_in[i, j] = 0.0  # Handle bins with zero total counts



            
            #### For region notIN

            Nuubar_highq_notIN = Nuubar_highq - Nuubar_highq_in
            Nuubar_highqbar_notIN = Nuubar_highqbar - Nuubar_highqbar_in
            Nddbar_highq_notIN = Nddbar_highq - Nddbar_highq_in
            Nddbar_highqbar_notIN = Nddbar_highqbar - Nddbar_highqbar_in

            NEvents_notIN = NEvents - NEvents_in
            N_non_udQ_notIN = N_non_udQ - N_non_udQ_in

            Nuubar_notIN = Nuubar_highq_notIN + Nuubar_highqbar_notIN
            Nddbar_notIN = Nddbar_highq_notIN + Nddbar_highqbar_notIN

            Nq_notIN = Nuubar_notIN + Nddbar_notIN + N_non_udQ_notIN

            if (Nuubar_highq_notIN + Nuubar_highqbar_notIN) > 0:
                Du_notIN[i, j] = (Nuubar_highq_notIN - Nuubar_highqbar_notIN)/(Nuubar_highq_notIN + Nuubar_highqbar_notIN)
            else:
                Du_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts

            if (Nddbar_highq_notIN + Nddbar_highqbar_notIN) > 0:
                Dd_notIN[i, j] = (Nddbar_highq_notIN - Nddbar_highqbar_notIN)/(Nddbar_highq_notIN + Nddbar_highqbar_notIN)
            else:
                Dd_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts

            if Nq_notIN > 0:
                Fu_Nq_notIN[i, j] = Nuubar_notIN/Nq_notIN
            else:
                Fu_Nq_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts

            if Nq_notIN > 0:
                Fd_Nq_notIN[i, j] = Nddbar_notIN/Nq_notIN
            else:
                Fd_Nq_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts

            if NEvents_notIN > 0:
                Fu_NEvents_notIN[i, j] = Nuubar_notIN/NEvents_notIN
            else:
                Fu_NEvents_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts

            if NEvents_notIN > 0:
                Fd_NEvents_notIN[i, j] = Nddbar_notIN/NEvents_notIN
            else:
                Fd_NEvents_notIN[i, j] = 0.0  # Handle bnotINs with zero total counts


            #### For region notOUT
            Nuubar_highq_notOUT = Nuubar_highq - Nuubar_highq_out
            Nuubar_highqbar_notOUT = Nuubar_highqbar - Nuubar_highqbar_out
            Nddbar_highq_notOUT = Nddbar_highq - Nddbar_highq_out
            Nddbar_highqbar_notOUT = Nddbar_highqbar - Nddbar_highqbar_out

            NEvents_notOUT = NEvents - NEvents_out
            N_non_udQ_notOUT = N_non_udQ - N_non_udQ_out

            Nuubar_notOUT = Nuubar_highq_notOUT + Nuubar_highqbar_notOUT
            Nddbar_notOUT = Nddbar_highq_notOUT + Nddbar_highqbar_notOUT

            Nq_notOUT = Nuubar_notOUT + Nddbar_notOUT + N_non_udQ_notOUT

            if (Nuubar_highq_notOUT + Nuubar_highqbar_notOUT) > 0:
                Du_notOUT[i, j] = (Nuubar_highq_notOUT - Nuubar_highqbar_notOUT)/(Nuubar_highq_notOUT + Nuubar_highqbar_notOUT)
            else:
                Du_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

            if (Nddbar_highq_notOUT + Nddbar_highqbar_notOUT) > 0:
                Dd_notOUT[i, j] = (Nddbar_highq_notOUT - Nddbar_highqbar_notOUT)/(Nddbar_highq_notOUT + Nddbar_highqbar_notOUT)
            else:
                Dd_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

            if Nq_notOUT > 0:
                Fu_Nq_notOUT[i, j] = Nuubar_notOUT/Nq_notOUT
            else:
                Fu_Nq_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

            if Nq_notOUT > 0:
                Fd_Nq_notOUT[i, j] = Nddbar_notOUT/Nq_notOUT
            else:
                Fd_Nq_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

            if NEvents_notOUT > 0:
                Fu_NEvents_notOUT[i, j] = Nuubar_notOUT/NEvents_notOUT
            else:
                Fu_NEvents_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

            if NEvents_notOUT > 0:
                Fd_NEvents_notOUT[i, j] = Nddbar_notOUT/NEvents_notOUT
            else:
                Fd_NEvents_notOUT[i, j] = 0.0  # Handle bnotOUTs with zero total counts

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

    np.save(f'{outputDir}/Du_in_{era}_{timeStamp}.npy', Du_in)
    np.save(f'{outputDir}/Dd_in_{era}_{timeStamp}.npy', Dd_in)
    np.save(f'{outputDir}/Fu_Nq_in_{era}_{timeStamp}.npy', Fu_Nq_in)
    np.save(f'{outputDir}/Fd_Nq_in_{era}_{timeStamp}.npy', Fd_Nq_in)
    np.save(f'{outputDir}/Fu_NEvents_in_{era}_{timeStamp}.npy', Fu_NEvents_in)
    np.save(f'{outputDir}/Fd_NEvents_in_{era}_{timeStamp}.npy', Fd_NEvents_in)

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
    np.save(f'{outputDir}/beta_ttz_edges_{era}_{timeStamp}.npy', beta_ttz_edges)


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

    Du_in_flat = Du_in.flatten()
    Dd_in_flat = Dd_in.flatten()
    Fu_Nq_in_flat = Fu_Nq_in.flatten()
    Fd_Nq_in_flat = Fd_Nq_in.flatten()
    Fu_NEvents_in_flat = Fu_NEvents_in.flatten()
    Fd_NEvents_in_flat = Fd_NEvents_in.flatten()

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
        'Du_in': Du_in_flat,
        'Dd_in': Dd_in_flat,
        'Fu_Nq_in': Fu_Nq_in_flat,
        'Fd_Nq_in': Fd_Nq_in_flat,
        'Fu_NEvents_in': Fu_NEvents_in_flat,
        'Fd_NEvents_in': Fd_NEvents_in_flat,
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