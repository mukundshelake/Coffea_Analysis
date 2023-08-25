import numpy as np, awkward as ak
import os
from coffea import processor
from coffea.analysis_tools import Weights
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema, NanoAODSchema
from coffea.util import load, save
import hist
from lib.helpers import getfileset


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
        events.Electron = events.Electron[
            (events.Electron.pt > 25) & (abs(events.Electron.eta) <= 3.0)
        ]
        req_ele = ak.count(events.Electron.pt, axis=1) >= 1
        events.Electron = ak.pad_none(events.Electron, 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 3.0)]
        req_jet = ak.count(events.Jet.pt, axis=1) >= 3
        events.Jet = ak.pad_none(events.Jet, 1)

        ## Other cuts
        ###### Add additional cuts here
        
        
        event_level = ak.fill_none(
            req_jet & req_ele, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        sele = events.Electron[event_level]
        sele = sele[:, 0]
        sjets = events.Jet[event_level]
        sjets = sjets[:, :3]
        
        print("No. of electrons: ", sum(event_level))
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
        output["electron_pt"].fill(ak.flatten(sele.pt, axis=-1), weight=weights.weight())
        output["electron_eta"].fill(ak.flatten(sele.eta, axis=-1), weight=weights.weight())      
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

   

inputDir = '/nfs/home/common/RUN2_UL/Tree_crab/SIXTEEN_preVFP'
# If you are using getFiles from the lib.helpers 
#fileset = {
#    "Data": getFiles('/home/mukund2/Projects/CoffeaTrial/Files/Data'),
#    "MC": getFiles('/home/mukund2/Projects/CoffeaTrial/Files/MC')
#}


# Generate the fileset in the appropriate format for the input directory.
fileset = getfileset(inputDir)        

iterative_run = processor.Runner(
    executor=processor.IterativeExecutor(compression=None),
    schema=NanoAODSchema,
    chunksize=200,
    maxchunks=10,
)

out = iterative_run(
    fileset,
    treename="Events",
    processor_instance=NanoProcessor(),
)
save(out, "outTest.coffea")  # save dictionary into coffea file
