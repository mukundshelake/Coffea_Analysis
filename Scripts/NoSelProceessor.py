import numpy as np
import hist
from coffea.analysis_tools import Weights
from coffea import processor
import awkward as ak
from coffea.nanoevents import NanoAODSchema
from coffea.lookup_tools import extractor
from coffea.analysis_tools import PackedSelection
import warnings
warnings.filterwarnings("ignore")


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
        ## create the dictionary contains 
        output = {
            "sumw" : processor.defaultdict_accumulator(float),
            "selEvents" : processor.defaultdict_accumulator(float) 
        }
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]

        ####################
        #    Selections    #
        ####################

        # selection = PackedSelection()

        event_level = np.full(len(events), True)
        if sum(event_level) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        #output["numEvents"] = sum(event_level)
        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(sum(event_level), storeIndividual=True)
        if isRealData:
            try:
                weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
            except:
                print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
        else:
            try:
                weights.add("PUWt", weight = events[event_level].puWeight, weightUp = events[event_level].puWeightUp, weightDown = events[event_level].puWeightDown)
            except:
                print(f"puWeight is not there for dataset {dataset}; adding 1s as puWeights")
                weights.add("PUWt", weight = np.ones(sum(event_level), dtype = float))
            
            try:
                weights.add("L1preFireWt", weight=events[event_level].L1PreFiringWeight.Nom, weightUp = events[event_level].L1PreFiringWeight.Up, weightDown = events[event_level].L1PreFiringWeight.Dn)
            except:
                print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                weights.add("L1preFireW", weight = np.ones(sum(event_level), dtype = float))
                
            try:
                weights.add("LHEWeightSign", weight = events[event_level].LHEWeight.originalXWGTUP/abs(events[event_level].LHEWeight.originalXWGTUP))
            except:
                print(f"LHEWeight is not there for dataset {dataset}; adding +1s as LHEWeightSign")
                weights.add("LHEWeightSign", weight = np.ones(sum(event_level), dtype = float))

        ####################
        #  Fill histogram  #
        ####################
        ## Filling the output dictionary
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.Generator.weight)
        output["selEvents"] = float(sum(event_level))

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator