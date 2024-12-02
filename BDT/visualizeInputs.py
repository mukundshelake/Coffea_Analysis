import pandas as pd
import matplotlib.pyplot as plt
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
    parser = argparse.ArgumentParser(description="Visualize input data from parquet files.")
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

    logging.info(f"Number of events: {len(df)}")

    # plot the ttbarpz distribution when y value is 0
    if 'ttbarpz' in df.columns:
        plt.hist(df[df['y'] == 2]['ttbarpz'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['ttbarpz'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['ttbarpz'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('ttbarpz')
        plt.xlim(-3000, 3000)
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/ttbarpz.png")
        logging.info(f"Plot saved to {outputDir}/ttbarpz.png")
        plt.clf()

    # plot similar for 'ttbar_mass'
    if 'ttbar_mass' in df.columns:
        plt.hist(df[df['y'] == 2]['ttbar_mass'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['ttbar_mass'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['ttbar_mass'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('ttbar_mass')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/ttbar_mass.png")
        logging.info(f"Plot saved to {outputDir}/ttbar_mass.png")
        plt.clf()

    # plot similar for FW1 binned between 2 * 10^6 to 1 * 10^6 with 100 bins
    if 'FW1' in df.columns:
        plt.hist(df[df['y'] == 2]['FW1'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['FW1'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['FW1'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('FW1')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/FW1.png")
        logging.info(f"Plot saved to {outputDir}/FW1.png")
        plt.clf()

    if 'FW2' in df.columns:
        plt.hist(df[df['y'] == 2]['FW2'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['FW2'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['FW2'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('FW2')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/FW2.png")
        logging.info(f"Plot saved to {outputDir}/FW2.png")
        plt.clf()

    if 'FW3' in df.columns:
        plt.hist(df[df['y'] == 2]['FW3'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['FW3'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['FW3'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('FW3')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/FW3.png")
        logging.info(f"Plot saved to {outputDir}/FW3.png")
        plt.clf()

    if 'FW4' in df.columns:
        plt.hist(df[df['y'] == 2]['FW4'], bins=100, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['FW4'], bins=100, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['FW4'], bins=100, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('FW4')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/FW4.png")
        logging.info(f"Plot saved to {outputDir}/FW4.png")
        plt.clf()

    if 'nJet' in df.columns:
        plt.hist(df[df['y'] == 2]['nJet'], bins=10, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['nJet'], bins=10, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['nJet'], bins=10, color = 'red', label='qq', density=True, histtype='step')
        plt.xlim(3, 12)
        plt.xlabel('nJet')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/nJet.png")
        logging.info(f"Plot saved to {outputDir}/nJet.png")
        plt.clf()

    if 'Jet_HT' in df.columns:
        plt.hist(df[df['y'] == 2]['Jet_HT'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Jet_HT'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Jet_HT'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlim(0, 2000)
        plt.xlabel('Jet_HT')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Jet_HT.png")
        logging.info(f"Plot saved to {outputDir}/Jet_HT.png")
        plt.clf()

    if 'pT_sum' in df.columns:
        plt.hist(df[df['y'] == 2]['pT_sum'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['pT_sum'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['pT_sum'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('pT_sum')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/pT_sum.png")
        logging.info(f"Plot saved to {outputDir}/pT_sum.png")
        plt.clf()

    if 'Sxx' in df.columns:
        plt.hist(df[df['y'] == 2]['Sxx'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Sxx'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Sxx'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Sxx')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Sxx.png")
        logging.info(f"Plot saved to {outputDir}/Sxx.png")
        plt.clf()

    if 'Syy' in df.columns:
        plt.hist(df[df['y'] == 2]['Syy'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Syy'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Syy'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Syy')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Syy.png")
        logging.info(f"Plot saved to {outputDir}/Syy.png")
        plt.clf()

    if 'Szz' in df.columns:
        plt.hist(df[df['y'] == 2]['Szz'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Szz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Szz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Szz')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Szz.png")
        logging.info(f"Plot saved to {outputDir}/Szz.png")
        plt.clf()

    if 'Sxy' in df.columns:
        plt.hist(df[df['y'] == 2]['Sxy'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Sxy'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Sxy'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Sxy')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Sxy.png")
        logging.info(f"Plot saved to {outputDir}/Sxy.png")
        plt.clf()

    if 'Sxz' in df.columns:
        plt.hist(df[df['y'] == 2]['Sxz'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Sxz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Sxz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Sxz')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Sxz.png")
        logging.info(f"Plot saved to {outputDir}/Sxz.png")
        plt.clf()

    if 'Syz' in df.columns:
        plt.hist(df[df['y'] == 2]['Syz'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['Syz'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['Syz'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('Syz')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/Syz.png")
        logging.info(f"Plot saved to {outputDir}/Syz.png")
        plt.clf()

    if 'S' in df.columns:
        plt.hist(df[df['y'] == 2]['S'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['S'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['S'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('S')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/S.png")
        logging.info(f"Plot saved to {outputDir}/S.png")
        plt.clf()

    if 'A' in df.columns:
        plt.hist(df[df['y'] == 2]['A'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['A'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['A'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('A')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/A.png")
        logging.info(f"Plot saved to {outputDir}/A.png")
        plt.clf()

    if 'AL' in df.columns:
        plt.hist(df[df['y'] == 2]['AL'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['AL'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['AL'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('AL')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/AL.png")
        logging.info(f"Plot saved to {outputDir}/AL.png")
        plt.clf()

    if 'P' in df.columns:
        plt.hist(df[df['y'] == 2]['P'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['P'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['P'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('P')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/P.png")
        logging.info(f"Plot saved to {outputDir}/P.png")
        plt.clf()

    if 'p2in' in df.columns:
        plt.hist(df[df['y'] == 2]['p2in'], bins=20, color='blue', label='qg', density=True, histtype='step')
        plt.hist(df[df['y'] == 1]['p2in'], bins=20, color = 'black', label='gg', density=True, histtype='step')
        plt.hist(df[df['y'] == 0]['p2in'], bins=20, color = 'red', label='qq', density=True, histtype='step')
        plt.xlabel('p2in')
        plt.ylabel('Number of events')
        plt.legend()
        plt.show()
        # save the figure
        plt.savefig(f"{outputDir}/p2in.png")
        logging.info(f"Plot saved to {outputDir}/p2in.png")
        plt.clf()

if __name__ == '__main__':
    main()
