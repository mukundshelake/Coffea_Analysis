import numpy as np, awkward as ak
import os, time
from coffea import processor
from coffea.analysis_tools import Weights
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema, NanoAODSchema
from coffea.util import load, save
import hist
from lib.helpers import getfileset, getFiles


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
            "muon_pt": hist.Hist(
                hist.axis.Regular(50, 0, 300, name="pt", label="$p_T$ [GeV]"),
                hist.storage.Weight(),
            ),
            "muon_eta": hist.Hist(
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
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Muon = events.Muon[
            (events.Muon.pt > 25) & (abs(events.Muon.eta) <= 3.0)
        ]
        req_mu = ak.count(events.Muon.pt, axis=1) >= 1
        events.Muon = ak.pad_none(events.Muon, 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 3.0)]
        req_jet = ak.count(events.Jet.pt, axis=1) >= 3
        events.Jet = ak.pad_none(events.Jet, 1)

        ## Other cuts
        ###### Add additional cuts here
	## HLT

        triggers = [
            "IsoMu24",
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
            req_jet & req_mu & req_trig, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        smu = events.Muon[event_level]
        smu = smu[:, 0]
        sjets = events.Jet[event_level]
        sjets = sjets[:, :3]
        
        print("No. of muons: ", sum(event_level))
        #output["numEvents"] = sum(event_level)
        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            genflavor = sjets.hadronFlavour
            weights.add("genweight", events[event_level].genWeight)
            genweiev = ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["pt"])[0])

        ####################
        #  Fill histogram  #
        ####################
        # The content filled to histogram should be flatten arrays(numpy arrays)
        output["muon_pt"].fill(ak.flatten(smu.pt, axis=-1), weight=weights.weight())
        output["muon_eta"].fill(ak.flatten(smu.eta, axis=-1), weight=weights.weight())      
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

   

era = 'EIGHTEEN'
lep = 'mu'
DataDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_{lep}'
MCDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/MC'

inputDir = [DataDir, MCDir]
# # If you are using getFiles from the lib.helpers 
# fileset = {
#    "Data": getfileset('/nfs/home/common/RUN2_UL/Tree_crab/SIXTEEN_preVFP/Data_el/'),
#    "MC": getfileset('/nfs/home/common/RUN2_UL/Tree_crab/SIXTEEN_preVFP/MC/')
# }


# Generate the fileset in the appropriate format for the input directory.
fileset = getfileset(inputDir)         

iterative_run = processor.Runner(
    executor=processor.IterativeExecutor(compression=None),
    schema=NanoAODSchema,
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
save(out, output_filename)  # save dictionary into coffea file
