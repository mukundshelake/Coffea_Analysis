#!/bin/bash

# Navigate to the NN folder
cd "$(dirname "$0")"

# Set the PYTHONPATH to include the parent directory
export PYTHONPATH=$(dirname "$PWD")

# Function to run a script and send a notification
run_script() {
    local script_name=$1
    local args=$2

    # Notify the start of the script
    python -c "from messege_bot import send_message; send_message('$script_name has started')"

    # Record the start time
    local start_time=$(date +%s)

    # Run the script
    python "$script_name" $args
    local status=$?

    # Record the end time
    local end_time=$(date +%s)

    # Calculate the total duration
    local duration=$((end_time - start_time))

    # Notify the end of the script with the total duration
    if [ $status -eq 0 ]; then
        python -c "from messege_bot import send_message; send_message('$script_name has completed successfully in $duration seconds')"
    else
        python -c "from messege_bot import send_message; send_message('$script_name encountered an error with status code $status in $duration seconds')"
    fi
}

# Check if a timestamp is provided as an argument
if [ -z "$1" ]; then
    # Generate a timestamp if not provided
    timestamp=$(date +"%Y%m%d_%H%M%S")
else
    timestamp=$1
fi

# Notify the start of the workflow
python -c "from messege_bot import send_message; send_message('Workflow started at $timestamp')"

# Record the start time of the workflow
workflow_start_time=$(date +%s)

# Run the classifier script
run_script "extractor.py" "-e UL2016preVFP -c ttbar_SemiLeptonic -t $timestamp -n 8"

# Uncomment and run additional scripts as needed
run_script "preprocesser.py" "-t $timestamp"
run_script "BDT.py" "-t $timestamp"
run_script "plotter.py" "-t $timestamp"
run_script "visualizeInputs.py" "-t $timestamp"

# Record the end time of the workflow
workflow_end_time=$(date +%s)

# Calculate the total duration of the workflow
workflow_duration=$((workflow_end_time - workflow_start_time))

# Notify the end of the workflow with the total duration
python -c "from messege_bot import send_message; send_message('Workflow completed in $workflow_duration seconds')"
