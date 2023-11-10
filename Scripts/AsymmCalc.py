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
            "noEvents":processor.defaultdict_accumulator(float),
            "t_eta": hist.Hist(
                        hist.axis.Regular(60, -3.0, 3.0, name="eta", label="$ \eta $ "),
                        hist.storage.Weight(),
                    ),
            "tbar_eta": hist.Hist(
                        hist.axis.Regular(60, -3.0, 3.0, name="eta", label="$ \eta $ "),
                        hist.storage.Weight(),
                    ),
        }
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]
        
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.Generator.weight)
        output["noEvents"] = len(events)

        if len(events) == 0:
            output["noEvents"] = 0.0
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        #output["numEvents"] = len(events)
        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events), storeIndividual=True)
        if isRealData:
            try:
                weights.add("L1preFireWt", weight=events.L1PreFiringWeight.Nom, weightUp = events.L1PreFiringWeight.Up, weightDown = events.L1PreFiringWeight.Dn)
            except:
                print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                weights.add("L1preFireW", weight = np.ones(len(events), dtype = float))
        else:
            try:
                weights.add("PUWt", weight = events.puWeight, weightUp = events.puWeightUp, weightDown = events.puWeightDown)
            except:
                print(f"puWeight is not there for dataset {dataset}; adding 1s as puWeights")
                weights.add("PUWt", weight = np.ones(len(events), dtype = float))
            
            try:
                weights.add("L1preFireWt", weight=events.L1PreFiringWeight.Nom, weightUp = events.L1PreFiringWeight.Up, weightDown = events.L1PreFiringWeight.Dn)
            except:
                print(f"L1preFireW is not there for dataset {dataset}; adding 1s as L1preFireW")
                weights.add("L1preFireW", weight = np.ones(len(events), dtype = float))
                
            try:
                weights.add("LHEWeightSign", weight = events.LHEWeight.originalXWGTUP/abs(events.LHEWeight.originalXWGTUP))
            except:
                print(f"LHEWeight is not there for dataset {dataset}; adding +1s as LHEWeightSign")
                weights.add("LHEWeightSign", weight = np.ones(len(events), dtype = float))
        ####################
        #  Fill histogram  #
        ####################
        ## Filling the output dictionary
        output["t_eta"].fill(ak.flatten(events.GenPart.eta[:,2], axis=-1), weight=weights.weight())
        output["t_eta"].fill(ak.flatten(events.GenPart.eta[:,3], axis=-1), weight=weights.weight())

        # print(output)
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator