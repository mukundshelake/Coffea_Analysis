import numpy as np, awkward as ak
import os, time
from coffea import processor
from coffea.analysis_tools import Weights
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema, NanoAODSchema
from coffea.util import load, save
import hist
from lib.helpers import getfileset, getFiles
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

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 35) & (abs(events.Electron.eta) <= 3.0)
        ]
        req_el = ak.count(events.Electron.pt, axis=1) >= 1
        events.Electron = ak.pad_none(events.Electron, 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 3.0)]
        req_jet = ak.count(events.Jet.pt, axis=1) >= 3
        events.Jet = ak.pad_none(events.Jet, 1)

        ## Other cuts
        ###### Add additional cuts here
	## HLT

        triggers = [
            "Ele32_eta2p1_WPTight_Gsf",
        ]
        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(np.array(triggers)[~checkHLT], " not exist in", dataset)
        trig_arrs = [
            events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)
        ]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        event_level = ak.fill_none(
            req_jet & req_el & req_trig, False)
        if len(events[event_level]) == 0:
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

        #### Initializing the extractor

        # ext = extractor()
        # ext.add_weight_sets(["* * /nfs/home/mukund/Projects/Coffea_Analysis/Scripts/SFs/UL2016_preVFP_Tight.root"])
        # ext.finalize()
        # evaluator = ext.make_evaluator()

        weights = Weights(len(events[event_level]), storeIndividual=True)

        if isRealData:
            print(f"working on data:{dataset}")
            genflavor = ak.zeros_like(sjets.pt)
        else:
            print(f"working on MC:{dataset}")
            genflavor = sjets.hadronFlavour
            if "QCD" not in dataset:
                weights.add("LHEweight", events[event_level].LHEWeight.originalXWGTUP/abs(events[event_level].LHEWeight.originalXWGTUP))
            genweiev = ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["pt"])[0])
            # eleSF = evaluator["EGamma_SF2D"](sel.eta, sel.pt)
            # eleSFerror = evaluator["EGamma_SF2D_error"](sel.eta, sel.pt)

            # weights.add(
            # "eleSF",
            # weight=ak.prod(eleSF),
            # weightUp=ak.prod(eleSF + eleSFerror),
            # )

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
# inputDir = [MCDir]

# Generate the fileset in the appropriate format for the input directory.
fileset = getfileset(inputDir)        

iterative_run = processor.Runner(
    executor=processor.IterativeExecutor(compression=None),
    schema=NanoAODSchema,
    # chunksize=5,
    # maxchunks=5,
)
out = iterative_run(
    fileset,
    treename="Events",
    processor_instance=NanoProcessor(),
    
)
# Get the base name of the script without the extension
script_filename = os.path.basename(__file__)
script_name = os.path.splitext(script_filename)[0]
# Generate a timestamp
timestamp = time.strftime("%Y%m%d_%H%M%S")

# Create a unique output name using script name and timestamp
output_filename = f"{script_name}_{timestamp}_{era}_{lep}.coffea"
outputfile = os.path.join(outputDir, output_filename)
save(out, outputfile)  # save dictionary into coffea file
