import json, os, argparse
import hist
import numpy as np
import dask
import awkward as ak
import hist.dask as hda
import dask_awkward as dak
from coffea.nanoevents import NanoAODSchema
from coffea import processor

from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
from coffea.util import save
import uproot
import logging

# Configure logger
logging.basicConfig(
    level=logging.INFO,  # Use INFO or WARNING in production
    format="%(asctime)s — %(levelname)s — %(message)s"
)
logger = logging.getLogger(__name__)

# Still suppress lower-level logs from noisy modules
for noisy_module in ["uproot", "dask", "fsspec", "urllib3"]:
    logging.getLogger(noisy_module).setLevel(logging.WARNING)

def remove_empty_files(fileset):
    cleaned_fileset = {}
    for dataset, content in fileset.items():
        files = content.get("files", {})
        valid_files = {}

        for filepath, treename in files.items():
            try:
                with uproot.open(f"{filepath}:{treename}") as tree:
                    if tree.num_entries > 0:
                        valid_files[filepath] = treename
                    else:
                        logger.warning(f"[Zero Entries] Skipping {filepath}")
            except Exception as e:
                logger.error(f"[Error] Skipping {filepath} due to: {e}")

        if valid_files:
            cleaned_fileset[dataset] = {"files": valid_files}
        else:
            logger.warning(f"[Warning] Skipping dataset '{dataset}' — all files empty or bad.")

    return cleaned_fileset


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        logger.info(f"Processing dataset: {dataset}")

        # Initialize histograms (fresh per dataset)
        logger.debug("Initializing histograms...")
        histograms = {
            "muon_pt": hist.Hist.new.Reg(40, 0.0, 400.0, name="pt", label="Leading Muon $p_T$ [GeV]").Double(),
            "muon_eta": hist.Hist.new.Reg(50, -2.5, 2.5, name="eta", label="Leading Muon $\eta$").Double(),
            "muon_phi": hist.Hist.new.Reg(64, -3.2, 3.2, name="phi", label="Leading Muon $\phi$").Double(),
            "muon_mass": hist.Hist.new.Reg(40, 0.0, 0.24, name="mass", label="Leading Muon Mass [GeV]").Double(),
            "jet_pt": hist.Hist.new.Reg(50, 0.0, 500.0, name="pt", label="Leading Jet $p_T$ [GeV]").Double(),
            "jet_eta": hist.Hist.new.Reg(50, -2.5, 2.5, name="eta", label="Leading Jet $\eta$").Double(),
            "jet_phi": hist.Hist.new.Reg(64, -3.2, 3.2, name="phi", label="Leading Jet $\phi$").Double(),
            "jet_mass": hist.Hist.new.Reg(50, 0.0, 100.0, name="mass", label="Leading Jet Mass [GeV]").Double(),
        }

        leading_muon_pt = events.Muon.pt[:, 0]
        leading_muon_eta = events.Muon.eta[:, 0]
        leading_muon_phi = events.Muon.phi[:, 0]
        leading_muon_mass = events.Muon.mass[:, 0]

        leading_jet_pt = events.Jet.pt[:,0]
        leading_jet_eta = events.Jet.eta[:,0]
        leading_jet_phi = events.Jet.phi[:,0]
        leading_jet_mass = events.Jet.mass[:,0]

        # --- Total Weights ---
        total_weight = ak.ones_like(leading_muon_pt)
        if hasattr(events, "MuonHLTWeight"):
            total_weight = total_weight *  events.MuonHLTWeight
        if hasattr(events, "MuonIDWeight"):
            total_weight = total_weight *  events.MuonIDWeight
        if hasattr(events, "LHEWeightSign"):
            total_weight = total_weight *  events.LHEWeightSign
        if hasattr(events, "bTaggingWeight"):
            total_weight = total_weight *  events.bTaggingWeight
        if hasattr(events, "L1PreFiringWeight_Nom"):
            total_weight = total_weight *  events.L1PreFiringWeight_Nom
        if hasattr(events, "puWeight"):
            total_weight = total_weight *  events.puWeight


        # Fill histograms
        histograms["muon_pt"].fill(pt=leading_muon_pt.compute(), weight=total_weight.compute())
        histograms["muon_eta"].fill(eta=leading_muon_eta.compute(), weight=total_weight.compute())
        histograms["muon_phi"].fill(phi=leading_muon_phi.compute(), weight=total_weight.compute())
        histograms["muon_mass"].fill(mass=leading_muon_mass.compute(), weight=total_weight.compute())
        histograms["jet_pt"].fill(pt=leading_jet_pt.compute(), weight=total_weight.compute())
        histograms["jet_eta"].fill(eta=leading_jet_eta.compute(), weight=total_weight.compute())
        histograms["jet_phi"].fill(phi=leading_jet_phi.compute(), weight=total_weight.compute())
        histograms["jet_mass"].fill(mass=leading_jet_mass.compute(), weight=total_weight.compute())

        return {
            "entries": ak.num(events, axis=0),
            "histos": histograms
        }


    def postprocess(self, accumulator):
        return accumulator


def main():
    parser = argparse.ArgumentParser(description="Process some eras.")
    allowed_eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']
    parser.add_argument('-e', '--era', choices=allowed_eras, required=True, help='Era to process')
    parser.add_argument('-s', '--sample', action='store_true')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag to include in output file name')
    args = parser.parse_args()

    logger.info(f"Selected era: {args.era}")
    logger.info(f"Sample mode: {args.sample}")
    logger.info(f"Output tag: {args.tag}")

    outputDir = "outputs"
    datasetFlag = 'sample' if args.sample else 'data'

    fileset = {}
    era = args.era
    with open(f'../../Datasets/selected_{datasetFlag}Files_{era}.json', 'r') as json_file:
        dicti = json.load(json_file)
        for pr in dicti['Data_mu']:
            datasetName = f'{era}_{pr}'
            fileset[datasetName] = {"files": dicti['Data_mu'][pr]}
        for pr in dicti['MC_mu']:
            datasetName = f'{era}_{pr}'
            fileset[datasetName] = {"files": dicti['MC_mu'][pr]}
    # print(f"Fileset: {fileset}")
    fileset = remove_empty_files(fileset)
    # print(f"Cleaned fileset: {fileset}")
    logger.info("Preprocessing fileset...")
    dataset_runnable, dataset_updated = preprocess(
        fileset,
        align_clusters=False,
        files_per_batch=1,
        skip_bad_files=True,
        save_form=False,
    )

    to_compute = apply_to_fileset(
        MyProcessor(),
        max_chunks(dataset_runnable, 300),
        schemaclass=NanoAODSchema,
    )

    (out,) = dask.compute(to_compute, scheduler='threads')

    # --- Save output .coffea file ---
    outputFile = f"{args.era}_{datasetFlag}_{args.tag}.coffea"
    save(out, os.path.join(outputDir, outputFile))
    logger.info(f"Saving output to {os.path.join(outputDir, outputFile)}")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()