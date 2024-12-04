import os
import dask
import numpy as np
import awkward as ak
import hist.dask as hda
from coffea import processor
from coffea.nanoevents import NanoEventsFactory,  BaseSchema, NanoAODSchema
from coffea.dataset_tools import apply_to_fileset, max_chunks, preprocess
import json, argparse
from coffea.util import save
from coffea.analysis_tools import PackedSelection


## Good Resource : https://indico.fnal.gov/event/11999/contributions/11335/attachments/7308/9405/JPilot_DPF2017.pdf


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        print(dataset)

        parts = dataset.split('_', 1)

        # Extract the era and channel
        era = parts[0]
        channel = parts[1]

        isData = False
        if 'Run' in channel:
            isData = True

        maps = {
            'btagThreshold': {
                'UL2016preVFP': 0.2598,
                'UL2016postVFP': 0.2489,
                'UL2017': 0.3040,
                'UL2018': 0.2783,
            }
        }

        selection = PackedSelection()
        # weight = Weights.weight()

        selection.add_multiple(
            {
                "atleastOneLep": ak.num(events.Muon) > 0,
                "atleastThreeJ": ak.num(events.Jet) > 2,
                "goodLeps" : ak.sum((events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId), axis=1) >= 1,
                "goodJets" : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4), axis = 1) >= 3,
                "BTag"  : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4) & (events.Jet.btagDeepFlavB > maps['btagThreshold'][era]), axis = 1) >= 2
            }
        )
        if 'UL2016' in dataset:
            selection.add("HLT", events.HLT.IsoTkMu24 | events.HLT.IsoMu24)
        elif 'UL2017' in dataset:
            selection.add("HLT", events.HLT.IsoMu27)
        else:
            selection.add("HLT", events.HLT.IsoMu24)
        # mask = selection.all("atleastOneLep", "atleastThreeJ")
        cutflow = selection.cutflow("atleastOneLep", "atleastThreeJ", "goodLeps", "goodJets", "BTag", "HLT")

        # honecut, hcutflow, labels = cutflow.yieldhist()

        results = cutflow.result()

        finalMask = results[-1][-1]

        nSelected = np.sum(finalMask.compute())
        print(f"nSelected: {nSelected}")

        sevents = events[finalMask]

        muons = events.Muon[(events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)][finalMask]
        muon_details = {
            "pt": muons.pt,
            "eta": muons.eta,
            "phi": muons.phi,
            "mass": muons.mass
        }
        # spt_flat = ak.flatten(muons.pt)
        # seta_flat = ak.flatten(muons.eta)
        # leading_muon = muons[ak.argmax(muons.pt, axis=1, keepdims=True)]
        # leading_muon_pt = ak.flatten(leading_muon.pt)
        # leading_muon_eta = ak.flatten(leading_muon.eta)
        # leading_muon_mass = ak.flatten(leading_muon.mass)
        # leading_muon_phi = ak.flatten(leading_muon.phi)
        # leading_muon_energy = np.sqrt(leading_muon_pt**2 * np.cosh(leading_muon_eta)**2 + leading_muon_mass**2)


        sjets = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4)][finalMask]
        jet_details = {
            "pt": sjets.pt,
            "eta": sjets.eta,
            "phi": sjets.phi,
            "mass": sjets.mass
        }
        # leadingJet = sjets[ak.argmax(sjets.pt, axis=1, keepdims=True)]
        # leadingJet_pt = ak.flatten(leadingJet.pt)
        # leadingJet_eta = ak.flatten(leadingJet.eta)
        # leadingJet_mass = ak.flatten(leadingJet.mass)
        # leadingJet_phi = ak.flatten(leadingJet.phi)
        # leadingJet_energy = np.sqrt(leadingJet_pt**2 * np.cosh(leadingJet_eta)**2 + leadingJet_mass**2)


        # subleadingJet = sjets[ak.argsort(sjets.pt, axis=1, ascending=False)[:, 1:2]]
        # subleadingJet_pt = ak.flatten(subleadingJet.pt)
        # subleadingJet_eta = ak.flatten(subleadingJet.eta)
        # subleadingJet_mass = ak.flatten(subleadingJet.mass)
        # subleadingJet_phi = ak.flatten(subleadingJet.phi)
        # subleadingJet_energy = np.sqrt(subleadingJet_pt**2 * np.cosh(subleadingJet_eta)**2 + subleadingJet_mass**2)

        # subsubleadingJet = sjets[ak.argsort(sjets.pt, axis=1, ascending=False)[:, 2:3]]
        # subsubleadingJet_pt = ak.flatten(subsubleadingJet.pt)
        # subsubleadingJet_eta = ak.flatten(subsubleadingJet.eta)
        # subsubleadingJet_mass = ak.flatten(subsubleadingJet.mass)
        # subsubleadingJet_phi = ak.flatten(subsubleadingJet.phi)
        # subsubleadingJet_energy = np.sqrt(subsubleadingJet_pt**2 * np.cosh(subsubleadingJet_eta)**2 + subsubleadingJet_mass**2)


        # Calculate px, py, pz for leading muon
        # leading_muon_px = leading_muon_pt * np.cos(leading_muon_phi)
        # leading_muon_py = leading_muon_pt * np.sin(leading_muon_phi)
        # leading_muon_pz = leading_muon_pt * np.sinh(leading_muon_eta)

        # Calculate px, py, pz for leading jet
        # leadingJet_px = leadingJet_pt * np.cos(leadingJet_phi)
        # leadingJet_py = leadingJet_pt * np.sin(leadingJet_phi)
        # leadingJet_pz = leadingJet_pt * np.sinh(leadingJet_eta)

        # Calculate px, py, pz for subleading jet
        # subleadingJet_px = subleadingJet_pt * np.cos(subleadingJet_phi)
        # subleadingJet_py = subleadingJet_pt * np.sin(subleadingJet_phi)
        # subleadingJet_pz = subleadingJet_pt * np.sinh(subleadingJet_eta)

        # Calculate px, py, pz for subsubleading jet
        # subsubleadingJet_px = subsubleadingJet_pt * np.cos(subsubleadingJet_phi)
        # subsubleadingJet_py = subsubleadingJet_pt * np.sin(subsubleadingJet_phi)
        # subsubleadingJet_pz = subsubleadingJet_pt * np.sinh(subsubleadingJet_eta)



        # Calculate elements of sphericity tensor
        # sxx = np.sum(leadingJet_px**2 + subleadingJet_px**2 + subsubleadingJet_px**2)
        # syy = np.sum(leadingJet_py**2 + subleadingJet_py**2 + subsubleadingJet_py**2)
        # szz = np.sum(leadingJet_pz**2 + subleadingJet_pz**2 + subsubleadingJet_pz**2)
        # sxy = np.sum(leadingJet_px*leadingJet_py + subleadingJet_px*subleadingJet_py + subsubleadingJet_px*subsubleadingJet_py)





        tops = ak.zip(
            {
                "pt" : sevents.GenPart.pt[:, 2],
                "eta": sevents.GenPart.eta[:, 2],
                "mass": sevents.GenPart.mass[:, 2],
                "phi" : sevents.GenPart.phi[:, 2]
            }
        )
        antitops = ak.zip(
            {
                "pt" : sevents.GenPart.pt[:, 3],
                "eta": sevents.GenPart.eta[:, 3],
                "mass": sevents.GenPart.mass[:, 3],
                "phi" : sevents.GenPart.phi[:, 3]
            }
        )
        
        p1_id = sevents.GenPart.pdgId[:, 0]
        p2_id = sevents.GenPart.pdgId[:, 1]

        # Get sevents 'nJet' by counting number of jets in jet_pt
        # nJet = ak.num(sjets.pt)

        # Calculate Jet_Ht 
        # Jet_HT = ak.sum(sjets.pt, axis=1)

        # calculate sum of Jet_HT leading muon pt and MET
        # pT_sum = Jet_HT + leading_muon_pt + sevents.MET.pt

        MET = sevents.MET.pt



        tpt = tops["pt"]
        teta = tops["eta"]
        # tmass = tops["mass"]
        # tphi = tops["phi"]
        tbarpt = antitops["pt"]
        tbareta = antitops["eta"]
        # tbarmass = antitops["mass"]
        # tbarphi = antitops["phi"]

        # Apply conditi

        # tpx = tpt*np.cos(tphi)
        # tpy = tpt*np.sin(tphi)
        tpz = tpt*np.sinh(teta)
        # tE = np.sqrt(tpt*tpt*np.cosh(teta)*np.cosh(teta) + tmass*tmass)

        # tbarpx = tbarpt*np.cos(tbarphi)
        # tbarpy = tbarpt*np.sin(tbarphi)
        tbarpz = tbarpt*np.sinh(tbareta)
        # tbarE = np.sqrt(tbarpt*tbarpt*np.cosh(tbareta)*np.cosh(tbareta) + tbarmass*tbarmass)

        # Four-momentum of the tbar system in lab frame
        # ttbarpx = tpx + tbarpx
        # ttbarpy = tpy + tbarpy
        ttbarpz = tpz + tbarpz
        # ttbarE = tE + tbarE

        # Boost velocity of the tbar system
        # beta_ttbar_x = ttbarpx/ttbarE
        # beta_ttbar_y = ttbarpy/ttbarE
        # beta_ttbar_z = ttbarpz/ttbarE

        # beta = np.sqrt(beta_ttbar_x*beta_ttbar_x + beta_ttbar_y*beta_ttbar_y + beta_ttbar_z*beta_ttbar_z)

        # gamma = 1.0 / np.sqrt(1 - beta**2)

        # bp = beta_ttbar_x*tpx + beta_ttbar_y*tpy + beta_ttbar_z*tpz
        # p_prime_x = tpx + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_x
        # p_prime_y = tpy + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_y
        # p_prime_z = tpz + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_z

        # Calculate the angle between the top quark and the z-axis in the tbar rest frame
        # cos_theta = p_prime_z/ np.sqrt(p_prime_x*p_prime_x + p_prime_y*p_prime_y + p_prime_z*p_prime_z)


        # Calculate the angle between the top quark and the z-axis in the tbar rest frame

        # mtt = np.sqrt(ttbarE*ttbarE - (ttbarpx*ttbarpx + ttbarpy*ttbarpy + ttbarpz*ttbarpz))

        # betattz = abs(ttbarpz)/ttbarE)
        params = {
            'p1_id': p1_id,
            'p2_id': p2_id,
            'ttbarpz': ttbarpz,
            # 'mtt': mtt,
            # 'cos_theta': cos_theta
            # 'leading_muon_px': leading_muon_px,
            # 'leading_muon_py': leading_muon_py,
            # 'leading_muon_pz': leading_muon_pz,
            # 'leading_muon_E': leading_muon_energy,
            # 'leadingJet_px': leadingJet_px,
            # 'leadingJet_py': leadingJet_py,
            # 'leadingJet_pz': leadingJet_pz,
            # 'leadingJet_E': leadingJet_energy,
            # 'subleadingJet_px': subleadingJet_px,
            # 'subleadingJet_py': subleadingJet_py,
            # 'subleadingJet_pz': subleadingJet_pz,
            # 'subleadingJet_E': subleadingJet_energy,
            # 'subsubleadingJet_px': subsubleadingJet_px,
            # 'subsubleadingJet_py': subsubleadingJet_py,
            # 'subsubleadingJet_pz': subsubleadingJet_pz,
            # 'subsubleadingJet_E': subsubleadingJet_energy,
            # 'nJet': nJet,
            # 'Jet_HT': Jet_HT,
            # 'pT_sum': pT_sum,
            'MET': MET,
            'jet_details': jet_details,
            'muon_details': muon_details
        }

        return {
            'nSelected': nSelected,
            'params': params
        }

    def postprocess(self, accumulator):
        pass


def split_fileset(fileset, num_chunks):
    files = fileset['UL2016preVFP_ttbar_SemiLeptonic']['files']
    items = list(files.items())
    mid = len(items) // 2
    chunk1 = dict(items[:mid])
    chunk2 = dict(items[mid:])

    files1 = {'UL2016preVFP_ttbar_SemiLeptonic': {'files': chunk1}}
    files2 = {'UL2016preVFP_ttbar_SemiLeptonic': {'files': chunk2}}
    return [files1, files2]


def main():
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

    # print(fileset)
    # exit()

    num_chunks = 2  # Adjust the number of chunks as needed
    fileset_chunks = split_fileset(fileset, num_chunks)

    for i, chunk in enumerate(fileset_chunks):
        dataset_runnable, dataset_updated = preprocess(
            chunk,
            align_clusters=False,
            files_per_batch=1,
            skip_bad_files=True,
            save_form=False,
        )

        to_compute = apply_to_fileset(
            MyProcessor(),
            max_chunks(dataset_runnable, 100),
            schemaclass=NanoAODSchema,
        )

        (out,) = dask.compute(to_compute, scheduler='threads')

        outputFile = f"BDTSkimmerOutput_{args.timestamp}_chunk{i}.coffea"
        save(out, os.path.join(outputDir, outputFile))
        print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

    # dataset_runnable, dataset_updated = preprocess(
    #     fileset,
    #     align_clusters=False,
    #     files_per_batch=1,
    #     skip_bad_files=True,
    #     save_form=False,
    # )

    # to_compute = apply_to_fileset(
    #     MyProcessor(),
    #     max_chunks(dataset_runnable, 100),
    #     schemaclass= NanoAODSchema,
    # )

    # (out,) = dask.compute(to_compute, scheduler='threads')
    

    # if not os.path.exists(outputDir):
    #     os.makedirs(outputDir)

    # outputFile = f"BDTSkimmerOutput_{args.timestamp}.coffea"
    # save(out, os.path.join(outputDir, outputFile))
    # print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
