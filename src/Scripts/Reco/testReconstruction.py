import json, os, argparse
import hist
import numpy as np
import dask
import awkward as ak
import hist.dask as hda
import dask_awkward as dak
from coffea.nanoevents import BaseSchema
from coffea import processor

from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
from coffea.util import save
import uproot
import logging
import yaml

# Configure logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s"
)
logger = logging.getLogger(__name__)

# Suppress lower-level logs from noisy modules
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


class RecoTestProcessor(processor.ProcessorABC):
    def __init__(self, era, config):
        self.era = era
        self.config = config
        # Load muon cuts from config
        self.muon_pt_lo = config['muon']['pt']['lo'][era]
        self.muon_abs_eta_hi = config['muon']['abs_eta']['hi'][era]
        self.muon_iso_hi = config['muon']['iso']['hi'][era]

    def process(self, events):
        dataset = events.metadata['dataset']
        logger.info(f"Processing dataset: {dataset}")

        # Filter events with successful reconstruction (chi2_status == 0)
        good_reco = events.chi2_status == 0
        # track number of bad reconstructions in this chunk
        total_before = int(ak.num(events, axis=0).compute())
        bad_reco_count = int(ak.sum(~good_reco).compute())
        events = events[good_reco]

        logger.info(f"Events: total={total_before}, good={total_before-bad_reco_count}, bad_reco={bad_reco_count}")

        # Initialize histograms
        logger.debug("Initializing histograms...")
        histograms = {
            "gen_top_mass": hist.Hist.new.Reg(60, 100.0, 400.0, name="mass", label="Gen Top Mass [GeV]").Double(),
            "gen_antitop_mass": hist.Hist.new.Reg(60, 100.0, 400.0, name="mass", label="Gen Anti-top Mass [GeV]").Double(),
            "gen_ttbar_mass": hist.Hist.new.Reg(70, 300.0, 1000.0, name="mass", label="Gen $t\\bar{t}$ Mass [GeV]").Double(),
            "reco_top_mass": hist.Hist.new.Reg(60, 100.0, 400.0, name="mass", label="Reco Top Mass [GeV]").Double(),
            "reco_antitop_mass": hist.Hist.new.Reg(60, 100.0, 400.0, name="mass", label="Reco Anti-top Mass [GeV]").Double(),
            "reco_ttbar_mass": hist.Hist.new.Reg(70, 300.0, 1000.0, name="mass", label="Reco $t\\bar{t}$ Mass [GeV]").Double(),
            "delta_top_mass": hist.Hist.new.Reg(60, -150.0, 150.0, name="mass", label="$\\Delta$ Top Mass (Gen - Reco) [GeV]").Double(),
            "delta_antitop_mass": hist.Hist.new.Reg(60, -150.0, 150.0, name="mass", label="$\\Delta$ Anti-top Mass (Gen - Reco) [GeV]").Double(),
            "delta_ttbar_mass": hist.Hist.new.Reg(80, -200.0, 200.0, name="mass", label="$\\Delta$ $t\\bar{t}$ Mass (Gen - Reco) [GeV]").Double(),
        }

        # --- Extract Gen-level masses ---
        logger.info("Extracting gen-level top quarks...")
        
        # Build GenPart 4-vectors from individual branches
        import vector
        vector.register_awkward()
        
        # Filter for last copy of top/antitop using statusFlags
        # isLastCopy flag is bit 13 (value 8192)
        isLastCopy_bit = 13
        has_lastcopy_flag = (events.GenPart_statusFlags & (1 << isLastCopy_bit)) != 0
        
        is_top = (events.GenPart_pdgId == 6) & has_lastcopy_flag
        is_antitop = (events.GenPart_pdgId == -6) & has_lastcopy_flag
        
        # Get top and antitop masses (first occurrence per event)
        gen_top_mass = ak.firsts(events.GenPart_mass[is_top])
        gen_antitop_mass = ak.firsts(events.GenPart_mass[is_antitop])
        
        # Debug: verify we're getting different top/antitop
        n_tops = int(ak.sum(ak.num(events.GenPart_mass[is_top], axis=1) > 0).compute())
        n_antitops = int(ak.sum(ak.num(events.GenPart_mass[is_antitop], axis=1) > 0).compute())
        logger.info(f"Gen particles: tops in {n_tops} events, antitops in {n_antitops} events")
        
        # Build 4-vectors for ttbar mass calculation
        gen_top_p4 = ak.zip({
            "pt": ak.firsts(events.GenPart_pt[is_top]),
            "eta": ak.firsts(events.GenPart_eta[is_top]),
            "phi": ak.firsts(events.GenPart_phi[is_top]),
            "mass": gen_top_mass
        }, with_name="Momentum4D")
        
        gen_antitop_p4 = ak.zip({
            "pt": ak.firsts(events.GenPart_pt[is_antitop]),
            "eta": ak.firsts(events.GenPart_eta[is_antitop]),
            "phi": ak.firsts(events.GenPart_phi[is_antitop]),
            "mass": gen_antitop_mass
        }, with_name="Momentum4D")
        
        # Calculate gen ttbar mass
        gen_ttbar_p4 = gen_top_p4 + gen_antitop_p4
        gen_ttbar_mass = gen_ttbar_p4.mass

        # --- Extract Reco-level masses with charge assignment ---
        logger.info("Extracting reco-level top quarks...")
        
        # Get leading muon with proper cuts from config
        muon_mask = (
            (events.Muon_pt > self.muon_pt_lo) & 
            (abs(events.Muon_eta) < self.muon_abs_eta_hi) & 
            events.Muon_tightId & 
            (events.Muon_pfRelIso04_all < self.muon_iso_hi)
        )
        
        # Get leading muon charge (first passing muon)
        muon_charge = ak.firsts(events.Muon_charge[muon_mask])
        
        # Build 4-vectors for Top_lep and Top_had
        top_lep_p4 = ak.zip({
            "pt": events.Top_lep_pt,
            "eta": events.Top_lep_eta,
            "phi": events.Top_lep_phi,
            "mass": events.Top_lep_mass
        }, with_name="Momentum4D")
        
        top_had_p4 = ak.zip({
            "pt": events.Top_had_pt,
            "eta": events.Top_had_eta,
            "phi": events.Top_had_phi,
            "mass": events.Top_had_mass
        }, with_name="Momentum4D")
        
        # Assign top/antitop based on muon charge
        # If muon charge > 0: leptonic side = top, hadronic side = antitop
        # If muon charge < 0: leptonic side = antitop, hadronic side = top
        is_positive = muon_charge > 0
        
        # Debug: check charge distribution
        n_positive = int(ak.sum(is_positive).compute())
        n_negative = int(ak.sum(~is_positive).compute())
        logger.info(f"Muon charges: positive={n_positive}, negative={n_negative}")
        
        reco_top_mass = ak.where(is_positive, top_lep_p4.mass, top_had_p4.mass)
        reco_antitop_mass = ak.where(is_positive, top_had_p4.mass, top_lep_p4.mass)
        
        # Reco ttbar mass
        reco_ttbar_p4 = top_lep_p4 + top_had_p4
        reco_ttbar_mass = reco_ttbar_p4.mass
        
        # Calculate delta masses (gen - reco)
        delta_top_mass = gen_top_mass - reco_top_mass
        delta_antitop_mass = gen_antitop_mass - reco_antitop_mass
        delta_ttbar_mass = gen_ttbar_mass - reco_ttbar_mass

        # --- Fill histograms ---
        logger.info("Filling histograms...")
        
        histograms["gen_top_mass"].fill(mass=gen_top_mass.compute())
        histograms["gen_antitop_mass"].fill(mass=gen_antitop_mass.compute())
        histograms["gen_ttbar_mass"].fill(mass=gen_ttbar_mass.compute())
        
        histograms["reco_top_mass"].fill(mass=reco_top_mass.compute())
        histograms["reco_antitop_mass"].fill(mass=reco_antitop_mass.compute())
        histograms["reco_ttbar_mass"].fill(mass=reco_ttbar_mass.compute())
        
        histograms["delta_top_mass"].fill(mass=delta_top_mass.compute())
        histograms["delta_antitop_mass"].fill(mass=delta_antitop_mass.compute())
        histograms["delta_ttbar_mass"].fill(mass=delta_ttbar_mass.compute())

        return {
            "entries": int(ak.num(events, axis=0).compute()),
            "histos": histograms,
            "bad_reco": bad_reco_count,
        }

    def postprocess(self, accumulator):
        return accumulator


def main():
    parser = argparse.ArgumentParser(description="Test ttbar reconstruction quality")
    allowed_eras = ['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018']
    parser.add_argument('-e', '--era', choices=allowed_eras, required=True, help='Era to process')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag to include in output file name')
    args = parser.parse_args()

    logger.info(f"Selected era: {args.era}")
    logger.info(f"Output tag: {args.tag}")

    # Load config file
    config_path = '/home/mukund/Projects/PhysicsTools/NanoAODTools/configs/masterConfig.yaml'
    logger.info(f"Loading config from {config_path}")
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    outputDir = "outputs"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # Load fileset for ttbar_SemiLeptonic only
    fileset = {}
    era = args.era
    json_path = f'/home/mukund/Projects/PhysicsTools/NanoAODTools/Datasets/{args.tag}_reco_{era}_dataFiles.json'
    
    logger.info(f"Loading fileset from {json_path}")
    with open(json_path, 'r') as json_file:
        dicti = json.load(json_file)
        if 'ttbar_SemiLeptonic' in dicti['MC_mu']:
            datasetName = f'{era}_ttbar_SemiLeptonic'
            fileset[datasetName] = {"files": dicti['MC_mu']['ttbar_SemiLeptonic']}
        else:
            logger.error("ttbar_SemiLeptonic not found in MC_mu!")
            return

    logger.info(f"Found {len(fileset)} dataset(s)")
    
    # Clean and preprocess fileset
    fileset = remove_empty_files(fileset)
    
    logger.info("Preprocessing fileset...")
    dataset_runnable, dataset_updated = preprocess(
        fileset,
        align_clusters=False,
        files_per_batch=1,
        skip_bad_files=True,
        save_form=False,
    )

    # Apply processor
    to_compute = apply_to_fileset(
        RecoTestProcessor(era=args.era, config=config),
        max_chunks(dataset_runnable, 300),
        schemaclass=BaseSchema,
    )

    logger.info("Computing results...")
    (out,) = dask.compute(to_compute, scheduler='threads')

    # Save output
    outputFile = f"{args.era}_ttbar_reco_test_{args.tag}.coffea"
    output_path = os.path.join(outputDir, outputFile)
    save(out, output_path)
    logger.info(f"Saved output to {output_path}")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
