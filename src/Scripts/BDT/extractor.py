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
import logging


## Good Resource : https://indico.fnal.gov/event/11999/contributions/11335/attachments/7308/9405/JPilot_DPF2017.pdf

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


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        logging.info("Initializing MyProcessor")
        pass

    def process(self, events):
        logging.info("Starting process method")
        dataset = events.metadata['dataset']
        logging.info(f"Processing dataset: {dataset}")

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
        logging.info(f"Number of selected events: {nSelected}")

        sevents = events[finalMask]

        muons = events.Muon[(events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.4) & (events.Muon.tightId)][finalMask]
        muon_details = {
            "pt": muons.pt,
            "eta": muons.eta,
            "phi": muons.phi,
            "mass": muons.mass
        }
        
        sjets = events.Jet[(events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.4)][finalMask]
        jet_details = {
            "pt": sjets.pt,
            "eta": sjets.eta,
            "phi": sjets.phi,
            "mass": sjets.mass
        }


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

        MET = sevents.MET.pt
        tpt = tops["pt"]
        teta = tops["eta"]
        tphi = tops["phi"]
        tmass = tops["mass"]
        tbarpt = antitops["pt"]
        tbareta = antitops["eta"]
        tbarphi = antitops["phi"]
        tbarmass = antitops["mass"]

        tpz = tpt*np.sinh(teta)
        tbarpz = tbarpt*np.sinh(tbareta)
        ttbarpz = tpz + tbarpz

        # Calculate t-tbar invariant mass
        t_energy = np.sqrt(tpt**2 * np.cosh(teta)**2 + tmass**2)
        tbar_energy = np.sqrt(tbarpt**2 * np.cosh(tbareta)**2 + tbarmass**2)
        ttbar_mass = np.sqrt((t_energy + tbar_energy)**2 - (tpt*np.cos(tphi) + tbarpt*np.cos(tbarphi))**2 - (tpt*np.sin(tphi) + tbarpt*np.sin(tbarphi))**2 - (tpz + tbarpz)**2)

        params = {
            'p1_id': p1_id,
            'p2_id': p2_id,
            'ttbarpz': ttbarpz,
            'MET': MET,
            'jet_details': jet_details,
            'muon_details': muon_details,
            'ttbar_mass': ttbar_mass
        }

        logging.info("Finished processing events")
        return {
            'nSelected': nSelected,
            'params': params
        }

    def postprocess(self, accumulator):
        logging.info("Postprocessing")
        pass


def split_fileset(fileset, num_chunks):
    logging.info("Splitting fileset")
    files = fileset['UL2016preVFP_ttbar_SemiLeptonic']['files']
    items = list(files.items())
    chunks = np.array_split(items, num_chunks)
    
    filesets = []
    for i, chunk in enumerate(chunks):
        chunk_dict = dict(chunk)
        filesets.append({'UL2016preVFP_ttbar_SemiLeptonic': {'files': chunk_dict}})
    
    return filesets


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

    parser.add_argument(
        '-n', '--num_chunks',
        type=int,
        default=2,
        help="Specify the number of chunks to split the fileset into."
    )

    args = parser.parse_args()

    outputDir = f"outputs/{args.timestamp}"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    setup_logging(os.path.splitext(os.path.basename(__file__))[0], outputDir)

    # Display the parsed arguments
    logging.info(f"Selected eras: {args.eras}")
    logging.info(f"Sample mode: {args.sample}")
    logging.info(f"Timestamp for book keeping: {args.timestamp}")
    logging.info(f"Number of chunks: {args.num_chunks}")
    if len(args.channels) > 0:
        logging.info(f"Channels: {args.channels}")
    else:
        logging.info("Channels: All")

    datasetFlag = 'data'

    if args.sample:
        datasetFlag = 'sample'

    fileset = {}
    for era in args.eras:
        with open(f'../../Datasets/{datasetFlag}Files_{era}.json', 'r') as json_file:
            dicti = json.load(json_file)
            for pr in dicti['Data_mu']:
                if args.onlyMC:
                    logging.info("Working on only the MC, ignoring data")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # logging.info(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['Data_mu'][pr]}
            for pr in dicti['MC_mu']:
                if args.onlyData:
                    logging.info("Working on only the data, ignoring MC")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # logging.info(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['MC_mu'][pr]}

    # logging.info(fileset)
    # exit()

    fileset_chunks = split_fileset(fileset, args.num_chunks)
    if args.sample:
        fileset_chunks = [fileset_chunks[0]]
    logging.info(fileset_chunks)

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
        if args.sample:
            outputFile = f"sample_BDTSkimmerOutput_{args.timestamp}_chunk{i}.coffea"
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        save(out, os.path.join(outputDir, outputFile))
        logging.info(f"Output file is stored in {os.path.join(outputDir, outputFile)}")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
