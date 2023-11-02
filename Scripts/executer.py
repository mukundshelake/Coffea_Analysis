import argparse
import os
import traceback
from datetime import datetime
from lib.helpers import getfileset, getFiles, logScript, scptoEOS
import warnings
warnings.filterwarnings("ignore")

from coffea.nanoevents import NanoAODSchema
from coffea import processor
from coffea.util import load, save

def executeAndLog(processor_address, meta_info=""):
    era = 'SIXTEEN_preVFP'
    lep = 'mu'
    DataDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/Data_{lep}'
    MCDir = f'/nfs/home/common/RUN2_UL/Tree_crab/{era}/MC'

    inputDir = [DataDir, MCDir] 

    # Generate the fileset in the appropriate format for the input directory.
    fileset = getfileset(inputDir) 
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
        iterative_run = processor.Runner(
            executor=processor.FuturesExecutor(compression=None, workers = 15),
            schema=NanoAODSchema,
            # chunksize=20,
            # maxchunks=10,
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
    parser.add_argument("-p", "--processor_address", required=True, help="Path to the processor to execute and log.")
    parser.add_argument("-m", "--meta_info", default="", help="Optional meta information.")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function with the specified script filename and meta information
    executeAndLog(args.processor_address, args.meta_info)
