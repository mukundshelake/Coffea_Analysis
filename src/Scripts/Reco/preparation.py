import os
import argparse
import uproot
import awkward as ak
import numpy as np
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
import h5py
import json


# Define physical constants
MW = 80.4  # W-boson mass (GeV)
MT = 172.5  # Top quark mass (GeV)
delta_MW = 2.085  # W-boson mass resolution (GeV)
delta_MT = 1.43  # Top quark mass resolution (GeV)

def create_output_path(input_file, output_directory):
    # Extract the last 4 subfolders from the input file path
    subfolders = os.path.normpath(input_file).split(os.sep)[-4:-1]
    # Create the new output file path
    output_file_path = os.path.join(output_directory, *subfolders, os.path.splitext(os.path.basename(input_file))[0] + '_output.h5')
    return output_file_path

def solve_npz(mpx, mpy, mpz, mpe, metx, mety):
    # Calculate the neutrino pz assuming a W-boson mass constraint
    alpha = MW**2 - mpe**2 + mpx**2 + mpy**2 + mpz**2 + 2*(mpx*metx + mpy*mety)
    a = mpz**2 - mpe**2
    b = alpha*mpz 
    c = alpha**2/4 - mpe**2*(metx**2 + mety**2)

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c
    # print(discriminant)
    # print(a)
    # print(b)
    # print(c)

    pz1 = ak.where(discriminant < 0, 0, (-b + np.sqrt(discriminant))/(2*a))
    pz2 = ak.where(discriminant < 0, 0, (-b - np.sqrt(discriminant))/(2*a))

    return pz1, pz2

def process_file(filePath, output_file):
    file = uproot.open(filePath)
    tree = file['Events']
    jets = ak.zip({
        'pt': tree['Jet_pt'].array(),
        'eta': tree['Jet_eta'].array(),
        'phi': tree['Jet_phi'].array(),
        'mass': tree['Jet_mass'].array(),
        'btag': tree['Jet_btagDeepFlavB'].array()
    })

    if len(jets.pt) < 1:
        print(f"No jets found in {filePath}. Skipping...")
        return


    sorted_jets = jets[ak.argsort(jets.btag, axis=1, ascending=False)]
    leading_bjets = sorted_jets[:, :2]
    non_bjets = sorted_jets[:, 2:]

    # Sort non_bjets by pt in descending order
    sorted_non_bjets = non_bjets[ak.argsort(non_bjets.pt, axis=1, ascending=False)]

    # Extract 2 leading non_bjets with highest pt
    leading_non_bjets = sorted_non_bjets[:, :2]

    bjet1 = leading_bjets[:, :1]
    bjet2 = leading_bjets[:, 1:2]
    ljet1 = leading_non_bjets[:, :1]
    ljet2 = leading_non_bjets[:, 1:2]

    bjet1_px = bjet1.pt * np.cos(bjet1.phi)
    bjet1_py = bjet1.pt * np.sin(bjet1.phi)
    bjet1_pz = bjet1.pt * np.sinh(bjet1.eta)
    bjet1_e = np.sqrt(bjet1_px**2 + bjet1_py**2 + bjet1_pz**2 + bjet1.mass**2)

    bjet2_px = bjet2.pt * np.cos(bjet2.phi)
    bjet2_py = bjet2.pt * np.sin(bjet2.phi)
    bjet2_pz = bjet2.pt * np.sinh(bjet2.eta)
    bjet2_e = np.sqrt(bjet2_px**2 + bjet2_py**2 + bjet2_pz**2 + bjet2.mass**2)

    ljet1_px = ljet1.pt * np.cos(ljet1.phi)
    ljet1_py = ljet1.pt * np.sin(ljet1.phi)
    ljet1_pz = ljet1.pt * np.sinh(ljet1.eta)
    ljet1_e = np.sqrt(ljet1_px**2 + ljet1_py**2 + ljet1_pz**2 + ljet1.mass**2)

    ljet2_px = ljet2.pt * np.cos(ljet2.phi)
    ljet2_py = ljet2.pt * np.sin(ljet2.phi)
    ljet2_pz = ljet2.pt * np.sinh(ljet2.eta)
    ljet2_e = np.sqrt(ljet2_px**2 + ljet2_py**2 + ljet2_pz**2 + ljet2.mass**2)

    jet_resolution = 0.05

    delta_bjet1_px = np.sqrt((jet_resolution*bjet1.pt * np.cos(bjet1.phi))**2)
    delta_bjet1_py = np.sqrt((jet_resolution*bjet1.pt * np.sin(bjet1.phi))**2)
    delta_bjet1_pz = np.sqrt((jet_resolution*bjet1.pt * np.sinh(bjet1.eta))**2)
    delta_bjet1_e = np.sqrt((bjet1.pt*bjet1.pt * jet_resolution*np.cosh(bjet1.eta)*np.cosh(bjet1.eta)/bjet1_e)**2)

    delta_bjet2_px = np.sqrt((jet_resolution*bjet2.pt * np.cos(bjet2.phi))**2)
    delta_bjet2_py = np.sqrt((jet_resolution*bjet2.pt * np.sin(bjet2.phi))**2)
    delta_bjet2_pz = np.sqrt((jet_resolution*bjet2.pt * np.sinh(bjet2.eta))**2)
    delta_bjet2_e = np.sqrt((bjet2.pt*bjet2.pt * jet_resolution*np.cosh(bjet2.eta)*np.cosh(bjet2.eta)/bjet2_e)**2)

    delta_ljet1_px = np.sqrt((jet_resolution*ljet1.pt * np.cos(ljet1.phi))**2)
    delta_ljet1_py = np.sqrt((jet_resolution*ljet1.pt * np.sin(ljet1.phi))**2)
    delta_ljet1_pz = np.sqrt((jet_resolution*ljet1.pt * np.sinh(ljet1.eta))**2)
    delta_ljet1_e = np.sqrt((ljet1.pt*ljet1.pt * jet_resolution*np.cosh(ljet1.eta)*np.cosh(ljet1.eta)/ljet1_e)**2)

    delta_ljet2_px = np.sqrt((jet_resolution*ljet2.pt * np.cos(ljet2.phi))**2)
    delta_ljet2_py = np.sqrt((jet_resolution*ljet2.pt * np.sin(ljet2.phi))**2)
    delta_ljet2_pz = np.sqrt((jet_resolution*ljet2.pt * np.sinh(ljet2.eta))**2)
    delta_ljet2_e = np.sqrt((ljet2.pt*ljet2.pt * jet_resolution*np.cosh(ljet2.eta)*np.cosh(ljet2.eta)/ljet2_e)**2)

    muons = ak.zip({
        'pt': tree['Muon_pt'].array(),
        'eta': tree['Muon_eta'].array(),
        'phi': tree['Muon_phi'].array(),
        'mass': tree['Muon_mass'].array(),
        'charge': tree['Muon_charge'].array(),
        'tightid': tree['Muon_tightId'].array(),
    })

    muon = muons[ak.argmax(muons.pt, axis=1, keepdims=True)]

    muon_px = muon.pt * np.cos(muon.phi)
    muon_py = muon.pt * np.sin(muon.phi)
    muon_pz = muon.pt * np.sinh(muon.eta)
    muon_e = np.sqrt(muon_px**2 + muon_py**2 + muon_pz**2 + muon.mass**2)
    muon_charge = muon.charge

    muon_resolution = 0.05

    delta_muon_px = np.sqrt((muon_resolution*muon.pt * np.cos(muon.phi))**2)
    delta_muon_py = np.sqrt((muon_resolution*muon.pt * np.sin(muon.phi))**2)
    delta_muon_pz = np.sqrt((muon_resolution*muon.pt * np.sinh(muon.eta))**2)
    delta_muon_e = np.sqrt((muon.pt*muon.pt * muon_resolution*np.cosh(muon.eta)*np.cosh(muon.eta)/muon_e)**2)

    mets = ak.zip({
        'pt': tree['MET_pt'].array(),
        'phi': tree['MET_phi'].array(),
    })

    met_px = mets.pt * np.cos(mets.phi)
    met_py = mets.pt * np.sin(mets.phi)


    met_pz1 , met_pz2 = solve_npz(muon_px, muon_py, muon_pz, muon_e, met_px, met_py)

    met_e1 = np.sqrt(met_px**2 + met_py**2 + met_pz1**2)
    met_e2 = np.sqrt(met_px**2 + met_py**2 + met_pz2**2)

    met_resolution = 0.05
    delta_met_px = np.sqrt((met_resolution*mets.pt * np.cos(mets.phi))**2)
    delta_met_py = np.sqrt((met_resolution*mets.pt * np.sin(mets.phi))**2)
    delta_met_pz1 = np.sqrt((met_resolution*mets.pt)**2)
    delta_met_pz2 = np.sqrt((met_resolution*mets.pt)**2)
    delta_met_e1 = np.sqrt((met_resolution*mets.pt)**2)
    delta_met_e2 = np.sqrt((met_resolution*mets.pt)**2)

    met_px_2d = ak.to_numpy(met_px).reshape(-1, 1)
    met_py_2d = ak.to_numpy(met_py).reshape(-1, 1)

    delta_met_e1_2d = ak.to_numpy(delta_met_e1).reshape(-1, 1)
    delta_met_e2_2d = ak.to_numpy(delta_met_e2).reshape(-1, 1)
    delta_met_px_2d = ak.to_numpy(delta_met_px).reshape(-1, 1)
    delta_met_py_2d = ak.to_numpy(delta_met_py).reshape(-1, 1)
    delta_met_pz1_2d = ak.to_numpy(delta_met_pz1).reshape(-1, 1)
    delta_met_pz2_2d = ak.to_numpy(delta_met_pz2).reshape(-1, 1)

    if 'SemiLeptonic' in filePath:
        top_pt = tree['GenPart_pt'].array()[:, 2]
        top_eta = tree['GenPart_eta'].array()[:, 2]
        top_phi = tree['GenPart_phi'].array()[:, 2]
        top_mass = tree['GenPart_mass'].array()[:, 2]

        top_px = top_pt * np.cos(top_phi)
        top_py = top_pt * np.sin(top_phi)
        top_pz = top_pt * np.sinh(top_eta)
        top_e = np.sqrt(top_px**2 + top_py**2 + top_pz**2 + top_mass**2)

        top_px_2d = ak.to_numpy(top_px).reshape(-1, 1)
        top_py_2d = ak.to_numpy(top_py).reshape(-1, 1)
        top_pz_2d = ak.to_numpy(top_pz).reshape(-1, 1)
        top_e_2d = ak.to_numpy(top_e).reshape(-1, 1)

        antitop_pt = tree['GenPart_pt'].array()[:, 3]
        antitop_eta = tree['GenPart_eta'].array()[:, 3]
        antitop_phi = tree['GenPart_phi'].array()[:, 3]
        antitop_mass = tree['GenPart_mass'].array()[:, 3]

        antitop_px = antitop_pt * np.cos(antitop_phi)
        antitop_py = antitop_pt * np.sin(antitop_phi)
        antitop_pz = antitop_pt * np.sinh(antitop_eta)
        antitop_e = np.sqrt(antitop_px**2 + antitop_py**2 + antitop_pz**2 + antitop_mass**2)

        antitop_px_2d = ak.to_numpy(antitop_px).reshape(-1, 1)
        antitop_py_2d = ak.to_numpy(antitop_py).reshape(-1, 1)
        antitop_pz_2d = ak.to_numpy(antitop_pz).reshape(-1, 1)
        antitop_e_2d = ak.to_numpy(antitop_e).reshape(-1, 1)

    muon_charge_2d = ak.to_numpy(muon_charge).reshape(-1, 1)



    perm1 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet1_e, bjet1_px, bjet1_py, bjet1_pz, 
                            bjet2_e, bjet2_px, bjet2_py, bjet2_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz, 
                            ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e1, met_px_2d, met_py_2d, met_pz1], axis=1)

    perm2 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet1_e, bjet1_px, bjet1_py, bjet1_pz,
                            bjet2_e, bjet2_px, bjet2_py, bjet2_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                            ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e2, met_px_2d, met_py_2d, met_pz2], axis=1)

    perm3 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet2_e, bjet2_px, bjet2_py, bjet2_pz,
                            bjet1_e, bjet1_px, bjet1_py, bjet1_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                            ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e1, met_px_2d, met_py_2d, met_pz1], axis=1)

    perm4 = np.concatenate([muon_e, muon_px, muon_py, muon_pz, bjet2_e, bjet2_px, bjet2_py, bjet2_pz,
                            bjet1_e, bjet1_px, bjet1_py, bjet1_pz, ljet1_e, ljet1_px, ljet1_py, ljet1_pz,
                            ljet2_e, ljet2_px, ljet2_py, ljet2_pz, met_e2, met_px_2d, met_py_2d, met_pz2], axis=1)

    res1 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz,
                            delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                            delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e1_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz1_2d], axis=1)

    res2 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz,
                            delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                            delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e2_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz2_2d], axis=1)

    res3 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz,
                            delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                            delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e1_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz1_2d], axis=1)

    res4 = np.concatenate([delta_muon_e, delta_muon_px, delta_muon_py, delta_muon_pz, delta_bjet2_e, delta_bjet2_px, delta_bjet2_py, delta_bjet2_pz,
                            delta_bjet1_e, delta_bjet1_px, delta_bjet1_py, delta_bjet1_pz, delta_ljet1_e, delta_ljet1_px, delta_ljet1_py, delta_ljet1_pz,
                            delta_ljet2_e, delta_ljet2_px, delta_ljet2_py, delta_ljet2_pz, delta_met_e2_2d, delta_met_px_2d, delta_met_py_2d, delta_met_pz2_2d], axis=1)
    
    if 'SemiLeptonic' in filePath:
        genInfo = np.concatenate([top_e_2d, top_px_2d, top_py_2d, top_pz_2d, antitop_e_2d, antitop_px_2d, antitop_py_2d, antitop_pz_2d, muon_charge_2d], axis=1)
    else:
        genInfo = np.concatenate([muon_charge_2d], axis=1)


    with h5py.File(output_file, 'w') as f:
        f.create_dataset('perm1', data=perm1)
        f.create_dataset('perm2', data=perm2)
        f.create_dataset('perm3', data=perm3)
        f.create_dataset('perm4', data=perm4)
        f.create_dataset('res1', data=res1)
        f.create_dataset('res2', data=res2)
        f.create_dataset('res3', data=res3)
        f.create_dataset('res4', data=res4)
        f.create_dataset('genInfo', data=genInfo)

    print(f"Saved permutations and resolutions to '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some eras.")

    # Define the allowed choices
    allowed_eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']

    # Add the --eras argument with choices and default value
    parser.add_argument(
        '-e', '--eras', 
        choices=allowed_eras, 
        nargs='*', 
        default=allowed_eras,
        help="Specify one or more eras. Allowed values are: 'UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018'. If not provided, all eras will be used by default."
    )

    parser.add_argument(
        '-c', '--channels',
        nargs='+',  # Accept one or more values
        type=str,
        default=[],  # Default to an empty list if no arguments are provided
        help="Specify one or more channels to process. By default, all channels"
    )


    # Add the --sample flag argument
    parser.add_argument(
        '-s', '--sample',
        action='store_true',
        help="If provided, the sample mode will be enabled."
    )

    parser.add_argument(
        '--skimmed',
        action='store_true',
        help="If provided, skimmed datasets will be used for analysis."
    )

    parser.add_argument(
        '--onlyData',
        action='store_true',
        help="If provided, only the data will be processed."
    )

    parser.add_argument(
        '--onlyMC',
        action='store_true',
        help="If provided, only the MC will be processed."
    )

    # Add the --output argument
    parser.add_argument(
        '-t','--timestamp',
        type=str,
        default='timestamp',
        help="Specify the timestamp'."
    )

    args = parser.parse_args()

    # Display the parsed arguments
    print(f"Selected eras: {args.eras}")
    print(f"Sample mode: {args.sample}")
    print(f"Timestamp for book keeping: {args.timestamp}")

    if len(args.channels) > 0:
        print(f"Channels: {args.channels}")
    else:
        print("Channels: All")


    outputDir = f"outputs/{args.timestamp}"
    datasetFlag = 'data'

    if args.sample:
        datasetFlag = 'sample'

    if args.skimmed:
        datasetFlag = 'skimmed_' + datasetFlag

    fileset = {}
    for era in args.eras:
        with open(f'../../Datasets/{datasetFlag}Files_{era}.json', 'r') as json_file:
            dicti = json.load(json_file)
            for pr in dicti['Data_mu']:
                if args.onlyMC:
                    print("Working on only the MC, ignoring data")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['Data_mu'][pr]}
            for pr in dicti['MC_mu']:
                if args.onlyData:
                    print("Working on only the data, ignoring MC")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['MC_mu'][pr]}

    output_directory = os.path.join('outputs', args.timestamp)
    os.makedirs(output_directory, exist_ok=True)

    input_files = []

    for dataset in fileset:
        input_files.extend(fileset[dataset]['files'])

    for input_file in input_files:
        print(f"Processing file: {input_file}")
        output_file = create_output_path(input_file, output_directory)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        process_file(input_file, output_file)


# Sample use: python preparation.py -t Jan20_2 -i /mnt/disk1/skimmed_Run2/UL2016preVFP/MC_mu/ttbar_SemiLeptonic/
