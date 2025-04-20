#!/usr/bin/env python3
"""
Workflow script that combines DataMCHist.py and histPlotter.py functionality
to process data and create plots in a single pipeline.
"""

import os
import argparse
import logging
import subprocess
import sys
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments for the workflow."""
    parser = argparse.ArgumentParser(
        description='Workflow for Coffea analysis and plotting',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Arguments from DataMCHist.py
    parser.add_argument('-e', '--era', 
                       choices=['UL2016preVFP', 'UL2016postVFP', 'UL2017', 'UL2018'],
                       required=True,
                       help='Era to process')
    parser.add_argument('-s', '--sample', 
                       action='store_true',
                       help='Run in sample mode')
    parser.add_argument('-t', '--tag', 
                       type=str,
                       required=True,
                       help='Tag to include in output file names')
    
    # Arguments from histPlotter.py
    parser.add_argument('--plots-dir',
                       default='plots',
                       help='Output directory for plots')
    parser.add_argument('--sample-info', 
                       default='sample_info.json',
                       help='Path to JSON file with cross sections and generated events')
    parser.add_argument('--skip-plotting',
                       action='store_true',
                       help='Skip the plotting step')
    
    return parser.parse_args()

def run_data_mc_hist(args):
    """Run the DataMCHist.py script with the given arguments."""
    logger.info("Running DataMCHist analysis...")
    
    # Build the command
    cmd = [
        sys.executable,  # Use the same Python interpreter
        str(Path(__file__).parent / 'DataMCHist.py'),
        '-e', args.era,
        '-t', args.tag
    ]
    
    if args.sample:
        cmd.append('-s')
    
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,  # Unbuffered
            universal_newlines=True
        )
        
        # Create file-like objects for stdout/stderr
        def reader(pipe, pipe_name):
            with pipe:
                for line in iter(pipe.readline, ''):
                    if pipe_name == 'stdout':
                        logger.info(line.rstrip())
                    else:
                        logger.warning(line.rstrip())
        
        # Start reader threads
        from threading import Thread
        stdout_thread = Thread(target=reader, args=(process.stdout, 'stdout'))
        stderr_thread = Thread(target=reader, args=(process.stderr, 'stderr'))
        stdout_thread.start()
        stderr_thread.start()
        
        # Wait for process to complete
        process.wait()
        
        # Wait for readers to finish
        stdout_thread.join()
        stderr_thread.join()
        
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)
            
    except subprocess.CalledProcessError as e:
        logger.error(f"DataMCHist.py failed with return code {e.returncode}")
        raise
    
    # Return the expected output file path
    dataset_flag = 'sample' if args.sample else 'data'
    output_file = f"{args.era}_{dataset_flag}_{args.tag}.coffea"
    return f'outputs/{output_file}'

def run_hist_plotter(args, input_file):
    """Run the histPlotter.py script with the given arguments."""
    if args.skip_plotting:
        logger.info("Skipping plotting step as requested")
        return
    
    logger.info("Running histPlotter...")
    
    # Create plots directory if it doesn't exist
    os.makedirs(args.plots_dir, exist_ok=True)
    
    # Build the command
    cmd = [
        sys.executable,  # Use the same Python interpreter
        str(Path(__file__).parent / 'histPlotter.py'),
        input_file,
        '--output-dir', args.plots_dir,
        '--sample-info', args.sample_info
    ]
    
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,  # Line buffered
            universal_newlines=True
        )
        
        # Stream stdout and stderr in real-time
        while True:
            output = process.stdout.readline() if process.stdout else ''
            error = process.stderr.readline() if process.stderr else ''
            
            if output == '' and error == '' and process.poll() is not None:
                break
                
            if output and output.strip():
                logger.info(output.strip())
            if error and error.strip():
                logger.warning(error.strip())
                
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)
            
    except subprocess.CalledProcessError as e:
        logger.error(f"histPlotter.py failed with return code {e.returncode}")
        raise

def main():
    args = parse_arguments()
    
    try:
        # Step 1: Run DataMCHist to process the data
        output_file = run_data_mc_hist(args)
        
        # Step 2: Run histPlotter to create plots
        if not args.skip_plotting:
            run_hist_plotter(args, output_file)
        
        logger.info("Workflow completed successfully")
    except Exception as e:
        logger.error(f"Workflow failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()