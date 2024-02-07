import argparse
import os, json
import traceback
from datetime import datetime
from lib.helpers import logScript, scptoEOS
import warnings
warnings.filterwarnings("ignore")

from coffea.nanoevents import NanoAODSchema
from coffea import processor
from coffea.util import load, save

def executeAndLog(processor_address, era, lep, istest, issample, isttbar, noLog, meta_info=""):
    # era = 'SIXTEEN_preVFP'
    # lep = 'el'
    # DataDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_{lep}'
    # MCDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/MC'

    # inputDir = [DataDir, MCDir] 

    # Generate the fileset in the appropriate format for the input directory.
    # fileset = getfileset(inputDir) 
    jsonPath = f'Datasets/dataFiles_{era}.json'

    if issample:
        jsonPath = f'Datasets/sampleFiles_{era}.json'

    with open(jsonPath, 'r') as json_file:
        fileLists = json.load(json_file) 
    
    fileset = {**fileLists[f'Data_{lep}'], **fileLists[f'MC_{lep}']}

    if isttbar:
        fileset = {"ttbar_SemiLeptonic": fileLists[f'MC_{lep}']["ttbar_SemiLeptonic"]}
    processorName = os.path.splitext(processor_address)[0].split('/')[-1]
    # Generate a timestamp
    startTime = datetime.now()
    timestamp = startTime.strftime('%b%d%Y_%H%M%S')

    # Generate unique coffea file, log file and EOSlog directory
    eosLogDir = startTime.strftime('%b%d%Y')
    coffeaOutput = f"../coffeaOutputs/{processorName}_{timestamp}_{era}_{lep}.coffea"
    logFile = f"../Logs/log_{processorName}_{timestamp}_{era}_{lep}.txt"
    errorFile = f"../Logs/errorLog_{processorName}_{timestamp}_{era}_{lep}.txt"
    runStatus = None
    print(processorName)
    try:
        exec(f"from {processorName} import NanoProcessor", globals())
        my_Processor = NanoProcessor()
        # Execute the script

        if istest:
            iterative_run = processor.Runner(
                executor=processor.FuturesExecutor(compression=None, workers = 15),
                schema=NanoAODSchema,
                chunksize= 20,
                maxchunks= 10,
            )
        else:
            iterative_run = processor.Runner(
                executor=processor.FuturesExecutor(compression=None, workers = 15),
                schema=NanoAODSchema
            )

        out = iterative_run(
            fileset,
            treename="Events",
            processor_instance=my_Processor, 
        )
        save(out, coffeaOutput) 
        runStatus = "Success"

    except Exception as e:
        print(e)
        traceback.print_exc()
        runStatus = "Failure"
        # Log the error to the error log file
        with open(errorFile, 'w') as error_log:
            error_log.write(f"Error message: {str(e)}\n")
            error_log.write("Stack trace:\n")
            traceback.print_exc(file=error_log)
    logScript(logFile, processor_address, coffeaOutput, meta_info, runStatus, errorFile)
    if not noLog:
        if runStatus=="Success":
            scptoEOS(coffeaOutput, eosLogDir)
        scptoEOS(logFile, eosLogDir)
        if runStatus=="Failure":
            scptoEOS(errorFile, eosLogDir)

# Specify the script file, output file, and optional meta information
if __name__ == "__main__":
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Log the processor content with error handling.")

    # Define command-line arguments for script filename and meta information
    parser.add_argument("-e", "--era", required=True, help="Data taking era to run on the analysisi on [SIXTEENpreVFP, SIXTEEN_postVFP, SEVENTEEN, EIGHTEEN]")
    parser.add_argument("-l", "--lep", required=True, help="Lepton to run analysis on")
    parser.add_argument("-p", "--processor_address", required=True, help="Path to the processor to execute and log.")
    parser.add_argument("-m", "--meta_info", default="", help="Optional meta information.")
    parser.add_argument('--test', action = "store_true", default= False, help="Is this a test run or full run")
    parser.add_argument('--sample', action = "store_true", default= False, help="Is this a sample run") 
    parser.add_argument('--ttbar', action = "store_true", default= False, help="Is this a only for ttbar MC samples")   
    parser.add_argument('--noLog', action = "store_true", default= False, help="Should we log the scripts to EOS")
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function with the specified script filename and meta information
    executeAndLog(args.processor_address, args.era, args.lep, args.test, args.sample, args.ttbar, args.noLog, args.meta_info)
