import numpy as np, awkward as ak
import os
from datetime import datetime
from coffea import processor
from coffea.analysis_tools import Weights
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema, NanoAODSchema
from coffea.util import load, save
import hist
from lib.helpers import getfileset, getFiles, logScript, scptoEOS
from coffea.lookup_tools import extractor


class NanoProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
        # Define histograms with hist class- https://hist.readthedocs.io/

    @property
    ## accumulating processed results
    def accumulator(self):
        return self._accumulator

    ## Place to process the Nanoevents, selections, weights, fill histogram is done here, considered as a loop with operating on columnar dataset
    def process(self, events):
        ## create the dictionary contains the
        output = {
            "sumw": processor.defaultdict_accumulator(float),
	        #"numEvents": processor.defaultdict_accumulator(float),
            "electron_pt": hist.Hist(
                hist.axis.Regular(50, 0, 300, name="pt", label="$p_T$ [GeV]"),
                hist.storage.Weight(),
            ),
            "electron_eta": hist.Hist(
                hist.axis.Regular(50, -3.0, 3.0, name="eta", label="$\eta$ "),
                hist.storage.Weight(),
            ),
            ## more than one dimension is possible
            "jet_pt": hist.Hist(
                hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour"),
                hist.axis.Regular(50, 0, 300, name="pt", label="$p_T$ [GeV]"),
                hist.storage.Weight(),
            ),
            "jet_eta": hist.Hist(
                hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour"),
                hist.axis.Regular(50, -3.0, 3.0, name="eta", label="$\eta$ [GeV]"),
                hist.storage.Weight(),
            ),             
        }
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.Generator.weight)

        ####################
        #    Selections    #
        ####################
        req_el = ak.count(events.Electron.pt, axis=1) >= 1
        req_jet = ak.count(events.Jet.pt, axis=1) >= 3
        event_level = ak.fill_none(req_jet & req_el, False)
        if len(events) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        sel = events.Electron[event_level]
        sel = sel[:, 0]
        sjets = events.Jet[event_level]
        sjets = sjets[:, :3]

        print("No. of electrons: ", sum(event_level))
        #output["numEvents"] = sum(event_level)
        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)

        if isRealData:
            print(f"working on data {dataset}")
            genflavor = ak.zeros_like(sjets.pt)
        else:
            print(f"working on MC {dataset}")
            genflavor = sjets.hadronFlavour
            weights.add("GeneratorWeight", events[event_level].Generator.weight)
            weights.add(
                "puWeight", 
                weight=events[event_level].puWeight,
                weightUp = events[event_level].puWeightUp,
                weightDown = events[event_level].puWeightDown)
            if "QCD" not in dataset:
                weights.add(
                    "alphaS",
                    weight = np.ones(len(events[event_level])),
                    # weightUp = events[event_level].LHEPdfWeight[:, 102],
                    # weightDown = events[event_level].LHEPdfWeight[:, 101],
                    )
            # weights.add("LHEScaleWeight",events[event_level].LHEScaleWeight)
            # weights.add("PSWeight",events[event_level].PSWeight)
            # weights.add("LHEPdfWeight",events[event_level].LHEPdfWeight)

        ####################
        #  Fill histogram  #
        ####################
        # The content filled to histogram should be flatten arrays(numpy arrays)
        output["electron_pt"].fill(ak.flatten(sel.pt, axis=-1), weight=weights.weight())
        output["electron_eta"].fill(ak.flatten(sel.eta, axis=-1), weight=weights.weight())      
        ## fill 2D histogram, since genweight is associate to the events, using ak.broadcast_arrays to assign the same weights for each jet
        output["jet_pt"].fill(
            ak.flatten(genflavor, axis=-1),
            ak.flatten(sjets.pt),
            weight=ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]),
        )
        output["jet_eta"].fill(
            ak.flatten(genflavor, axis=-1),
            ak.flatten(sjets.eta),
            weight=ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["eta"])[0]),
        )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator

   

era = 'SIXTEEN_preVFP'
lep = 'el'
DataDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_{lep}'
MCDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/MC'
outputDir = f'../coffeaOutputs'

inputDir = [DataDir, MCDir] 

# Generate the fileset in the appropriate format for the input directory.
fileset = getfileset(inputDir)        

iterative_run = processor.Runner(
    executor=processor.IterativeExecutor(compression=None, workers = 4),
    schema=NanoAODSchema,
    # chunksize = 20,
    # maxchunks = 10,
)

out = iterative_run(
    fileset,
    treename="Events",
    processor_instance=NanoProcessor(),
    
)
# Get the base name of the script without the extension
script_filename = os.path.abspath(__file__)
script_name = os.path.splitext(script_filename)[0].split('/')[-1]
# Get the current date and time
current_datetime = datetime.now()

# Format the date in 'AbbDDYYYY' format
timestamp = current_datetime.strftime('%b%d%Y_%H%M%S')
eosLogDir = current_datetime.strftime('%b%d%Y')
# Create a unique output name using script name and timestamp
output_filename = f"{script_name}_{timestamp}_{era}_{lep}.coffea"
outputfile = os.path.join(outputDir, output_filename)
save(out, outputfile)  # save dictionary into coffea file

# Log the stuff
metaTxt = "Run: Adding LHEpdfWeights weights and PU up/down variations to (generator weights + PU weights on 1 lep and 3 jets)"
logScript("logger.txt", script_filename, outputfile, metaTxt)
scptoEOS(outputfile, eosLogDir)
