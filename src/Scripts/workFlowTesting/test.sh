#!/bin/bash

# Check if running on sample data
if [ "$1" == "sample" ]; then
    run_sample=true
    sample_indicator="_sample"
    skimmer_args="-s"
else
    run_sample=false
    sample_indicator=""
    skimmer_args=""
fi

# Generate a timestamp for the log filename and output filename
timestamp=$(date +"%Y-%m-%d_%H-%M-%S${sample_indicator}")
log_file="script_log_${timestamp}.log"
output_file="output_${timestamp}.coffea"

# Function to log messages with an optional timestamp
log() {
    local message=$1
    local with_timestamp=${2:-true}  # Default to true if not specified

    if [ "$with_timestamp" = true ]; then
        local current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
        echo "$current_timestamp - $message" >> "$log_file"
    else
        echo "$message" >> "$log_file"
    fi
}

# Function to get the current Git commit hash
get_git_commit() {
    # Try to get the Git commit hash
    git_commit=$(git rev-parse HEAD 2>/dev/null)
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "$git_commit"
    else
        echo "No Git Repository"
    fi
}

# Add and commit files to Git
add_and_commit_to_git() {
    local commit_message=$1

    # Add files to Git (modify this to specify the files you want to add)
    git add .

    # Commit files to Git
    git commit -m "$commit_message"
}

# Example usage
log "Script started"

# Add and commit files to Git
add_and_commit_to_git "Automated commit before running script"
log "Files added and committed to Git" false  # No timestamp for this message

# Get the current Git commit hash
commit_hash=$(get_git_commit)
log "Current Git commit hash: $commit_hash" false  # No timestamp for this message

# Don't touch anything above
# -------------------------------------------------------------------------------------------
# Your script commands below
log "Running command: python newskimmer.py $skimmer_args -o $output_file"
python newskimmer.py $skimmer_args -o $output_file 2>&1 | tee -a "$log_file"

log "Running command: python messageBot.py"
python messegeBot.py 2>&1 | tee -a "$log_file"




#--------------------------------------------------------------------------------------------------
log "Script completed"
