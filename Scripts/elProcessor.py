import numpy as np, awkward as ak
import os, time
from coffea import processor
from coffea.analysis_tools import Weights
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema, NanoAODSchema
from coffea.util import load, save
import hist
from lib.helpers import getfileset, getFiles, logScript
from coffea.lookup_tools import extractor
import json
##############################################################################3
## Things to do
#1. implement as much config as possible
#2. Increase the cores
#3. Add config.json to write function to logger

#################################################################################
# Load the configuration from config.json
with open("configEL.json", "r") as config_file:
    config = json.load(config_file)

era = config["era"]
lep = config["lep"]
DataDir = config["DataDir"]
MCDir = config["MCDir"]
outputDir = config["outputDir"]
inputDir = config["inputDir"]

# Load numerical values from the configuration
electron_pt_threshold = config["electron_pt_threshold"]
jet_pt_threshold = config["jet_pt_threshold"]
jet_eta_threshold = config["jet_eta_threshold"]
chunksize = config["chunkSize"]
maxchunk = config["maxChunks"]

##################################################################################
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
            "sumw": processor.defaultdict_accumulator(float)           
        }
        isRealData = not hasattr(events, "Generator")
        dataset = events.metadata["dataset"]
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.Generator.weight)
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator

# Generate the fileset in the appropriate format for the input directory.
fileset = getfileset(inputDir)        

iterative_run = processor.Runner(
    executor=processor.IterativeExecutor(compression=None),
    schema=NanoAODSchema,
    chunksize= chunksize,
    maxchunks= maxchunk,
)

out = iterative_run(
    fileset,
    treename="Events",
    processor_instance=NanoProcessor(),
    
)
# Get the base name of the script without the extension
script_filename = os.path.abspath(__file__)
script_name = os.path.splitext(script_filename)[0]
# Generate a timestamp
timestamp = time.strftime("%Y%m%d_%H%M%S")

# Create a unique output name using script name and timestamp
output_filename = f"{script_name}_{timestamp}_{era}_{lep}.coffea"
outputfile = os.path.join(outputDir, output_filename)
save(out, outputfile)  # save dictionary into coffea file


# Log the stuff
metaTxt = "Config trial with chunkSize and maxchunk"
logScript("logger.txt", script_filename, outputfile, metaTxt)