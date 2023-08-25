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
            "noEvents": processor.defaultdict_accumulator(float)
        }
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
            output["noEvents"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)
            output["noEvents"] = len(events)

        ####################
        #    No Selections    #
        ####################
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator

inputDir = '/nfs/home/common/RUN2_UL/Tree_crab/SIXTEEN_preVFP'


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
save(out, "outNoSel.coffea")  # save dictionary into coffea file

from coffea.util import load
output = load("outNoSel.coffea")
print(output)
