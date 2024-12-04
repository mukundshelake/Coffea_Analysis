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
from coffea.analysis_tools import PackedSelection, Weights

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

        selection = PackedSelection()
        # weight = Weights.weight()

        maps = {
            'btagThreshold': {
                'UL2016preVFP': 0.2598,
                'UL2016postVFP': 0.2489,
                'UL2017': 0.3040,
                'UL2018': 0.2783,
            }
        }

        selection.add_multiple(
            {
                "atleastOneLep": ak.num(events.Muon) > 0,
                "atleastThreeJ": ak.num(events.Jet) > 2,
                "goodLeps" : ak.sum((events.Muon.pt >= 35.0) & (abs(events.Muon.eta) <= 2.1) & (events.Muon.tightId), axis=1) >= 1,
                "goodJets" : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.5), axis = 1) >= 3,
                "BTag"  : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.5) & (events.Jet.btagDeepFlavB > maps['btagThreshold'][era]), axis = 1) >= 3
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
        results = cutflow.result()

        finalMask = results[-1][-1]

        nSelected = np.sum(finalMask.compute())
        tops = ak.zip(
            {
                "pt" : events[finalMask].GenPart.pt[:, 2],
                "eta": events[finalMask].GenPart.eta[:, 2],
                "mass": events[finalMask].GenPart.mass[:, 2],
                "phi" : events[finalMask].GenPart.phi[:, 2]
            }
        )
        antitops = ak.zip(
            {
                "pt" : events[finalMask].GenPart.pt[:, 3],
                "eta": events[finalMask].GenPart.eta[:, 3],
                "mass": events[finalMask].GenPart.mass[:, 3],
                "phi" : events[finalMask].GenPart.phi[:, 3]
            }
        )
        
        p1_id = events[finalMask].GenPart.pdgId[:, 0]
        p2_id = events[finalMask].GenPart.pdgId[:, 1]

        p1ID = ak.where(p1_id == 21, 0, p1_id)
        p2ID = ak.where(p2_id == 21, 0, p2_id)

        tpt = tops["pt"]
        teta = tops["eta"]
        tmass = tops["mass"]
        tphi = tops["phi"]
        yt = abs(teta - 0.50*np.tanh(teta)*np.square(tmass/tpt))
        tbarpt = antitops["pt"]
        tbareta = antitops["eta"]
        tbarmass = antitops["mass"]
        tbarphi = antitops["phi"]
        ytbar = abs(tbareta - 0.50*np.tanh(tbareta)*np.square(tbarmass/tbarpt))

        # Apply conditi

        tpx = tpt*np.cos(tphi)
        tpy = tpt*np.sin(tphi)
        tpz = tpt*np.sinh(teta)
        tE = np.sqrt(tpt*tpt*np.cosh(teta)*np.cosh(teta) + tmass*tmass)

        deltay = yt - ytbar

        tbarpx = tbarpt*np.cos(tbarphi)
        tbarpy = tbarpt*np.sin(tbarphi)
        tbarpz = tbarpt*np.sinh(tbareta)
        tbarE = np.sqrt(tbarpt*tbarpt*np.cosh(tbareta)*np.cosh(tbareta) + tbarmass*tbarmass)

        # Four-momentum of the tbar system in lab frame
        ttbarpx = tpx + tbarpx
        ttbarpy = tpy + tbarpy
        ttbarpz = tpz + tbarpz
        ttbarE = tE + tbarE

        # Boost velocity of the tbar system
        beta_ttbar_x = ttbarpx/ttbarE
        beta_ttbar_y = ttbarpy/ttbarE
        beta_ttbar_z = ttbarpz/ttbarE

        beta = np.sqrt(beta_ttbar_x*beta_ttbar_x + beta_ttbar_y*beta_ttbar_y + beta_ttbar_z*beta_ttbar_z)

        gamma = 1.0 / np.sqrt(1 - beta**2)

        bp = beta_ttbar_x*tpx + beta_ttbar_y*tpy + beta_ttbar_z*tpz
        p_prime_x = tpx + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_x
        p_prime_y = tpy + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_y
        p_prime_z = tpz + ((gamma - 1) * bp / beta**2 - gamma * tE) * beta_ttbar_z

        # Calculate the angle between the top quark and the z-axis in the tbar rest frame
        cos_theta = p_prime_z/ np.sqrt(p_prime_x*p_prime_x + p_prime_y*p_prime_y + p_prime_z*p_prime_z)


        # Calculate the angle between the top quark and the z-axis in the tbar rest frame

        mtt = np.sqrt(ttbarE*ttbarE - (ttbarpx*ttbarpx + ttbarpy*ttbarpy + ttbarpz*ttbarpz))

        # betattz = abs(ttbarpz)/ttbarE
        
        yt2D = (
            hda.Hist.new
            .Reg(4, -1.0, 1.0, label = "$c*$", name = "c")
            .Reg(24, 300, 1500, label = "$m_tt$", name = "m_tt")
            # .Reg(4, 0, 1.00, label = "$beta_ttz$", name = "beta_ttz")
            .Reg(4, 0, 2.4, label="$y_t$", name = "y_t")
            .Reg(4, 0, 2.4, label="$y_tbar$", name = "y_tbar")
            .Reg(2, -2.4, 2.4, label='$deltay$', name='deltay')
            .Integer(-5, 6, label="p1", name="p1")
            .Integer(-5, 6, label="p2", name="p2")
            .Double()
        )
        yt2D.fill(
            c = cos_theta,
            m_tt = mtt,
            # beta_ttz = betattz,
            y_t = yt,
            y_tbar = ytbar,
            deltay = deltay,
            p1 = p1ID,
            p2 = p2ID
        )
        return {
                "entries": ak.num(events, axis=0),
                "yMatrix": yt2D,
            }

    def postprocess(self, accumulator):
        pass

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
            for pr in dicti['Data_el']:
                if args.onlyMC:
                    print("Working on only the MC, ignoring data")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['Data_el'][pr]}
            for pr in dicti['MC_el']:
                if args.onlyData:
                    print("Working on only the data, ignoring MC")
                    continue
                if len(args.channels) > 0:
                    if pr not in args.channels:
                        # print(f'skipping {era}_{pr} as not in list')
                        continue
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['MC_el'][pr]}

    # print(fileset)
    # exit()
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
        schemaclass= NanoAODSchema,
    )

    (out,) = dask.compute(to_compute, scheduler='threads')
    

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    outputFile = f"LHSskimmerOutput_{args.timestamp}.coffea"
    save(out, os.path.join(outputDir, outputFile))
    print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
