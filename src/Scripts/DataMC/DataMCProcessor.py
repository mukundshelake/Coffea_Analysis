import json, os, argparse
import hist
import dask
import awkward as ak
import hist.dask as hda
import dask_awkward as dak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
from distributed import Client
from coffea.analysis_tools import PackedSelection, Weights
from coffea.util import save


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def process(self, events):
        dataset = events.metadata['dataset']
        print(dataset)

        selection = PackedSelection()
        # weight = Weights.weight()

        selection.add_multiple(
            {
                "atleastOneLep": ak.num(events.Electron) > 0,
                "atleastThreeJ": ak.num(events.Jet) > 2,
                "goodLeps" : ak.sum((events.Electron.pt >= 35.0) & (abs(events.Electron.eta) <= 2.1) & (events.Electron.cutBased == 4), axis=1) >= 1,
                "goodJets" : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.5), axis = 1) >= 3,
                "BTag"  : ak.sum((events.Jet.pt >= 30.0) & (abs(events.Jet.eta) < 2.5) & (events.Jet.btagDeepFlavB > 0.2598), axis = 1) >= 3
            }
        )

        if 'UL2017' in dataset:
            selection.add("HLT", events.HLT.Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned | events.HLT.Ele35_WPTight_Gsf)
        elif 'UL2018' in dataset:
            selection.add("HLT", events.HLT.Ele32_WPTight_Gsf)
        else:
            selection.add("HLT", events.HLT.Ele32_eta2p1_WPTight_Gsf)


        # mask = selection.all("atleastOneLep", "atleastThreeJ")
        cutflow = selection.cutflow("atleastOneLep", "atleastThreeJ", "goodLeps", "goodJets", "BTag", "HLT")

        honecut, hcutflow, labels = cutflow.yieldhist()

        return {
            dataset: {
                "entries": ak.num(events, axis=0),
                "honecut" : honecut,
                "hcutflow": hcutflow,
                "label": labels
            }
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

    # Add the --sample flag argument
    parser.add_argument(
        '-s', '--sample',
        action='store_true',
        help="If provided, the sample mode will be enabled."
    )

    # Add the --output argument
    parser.add_argument(
        '-o','--output',
        type=str,
        default='output.coffea',
        help="Specify the output file name. Default is 'output.coffea'."
    )

    args = parser.parse_args()

    # Display the parsed arguments
    print(f"Selected eras: {args.eras}")
    print(f"Sample mode: {args.sample}")
    print(f"Output file: {args.output}")

    outputDir = "../outputs"
    datasetFlag = 'data'

    if args.sample:
        datasetFlag = 'sample'

    fileset = {}
    for era in args.eras:
        with open(f'../../Datasets/{datasetFlag}Files_{era}.json', 'r') as json_file:
            dicti = json.load(json_file)
            for pr in dicti['Data_el']:
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['Data_el'][pr]}
            for pr in dicti['MC_el']:
                datasetName = f'{era}_{pr}'
                fileset[datasetName] = {"files": dicti['MC_el'][pr]}

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
    
    outputFile = args.output
    save(out, os.path.join(outputDir, outputFile))
    print(f"Output file is stored in {os.path.join(outputDir, outputFile)}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()