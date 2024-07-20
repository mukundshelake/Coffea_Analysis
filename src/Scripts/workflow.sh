############################################################################
#!/bin/bash

# Initialize variables
run_sample=false
sample_indicator=""
skimmer_args=""
description=""

# Check if running on sample data
if [ "$1" == "sample" ]; then
    run_sample=true
    sample_indicator="_sample"
    skimmer_args="-s"
fi

# Check if a description is provided
if [ -n "$2" ]; then
    description=$2
fi

# Generate a timestamp for the log filename and output filename
timestamp=$(date +"%Y%m%d_%H%M%S${sample_indicator}")
log_file="log_${timestamp}.log"
output_file="output_${timestamp}.coffea"

# Function to log messages with an optional timestamp
log() {
    local message=$1
    local with_timestamp=${2:-true}  # Default to true if not specified

    if [ "$with_timestamp" = true ]; then
        local current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
        echo "$current_timestamp - $message" | tee -a "$log_file"
    else
        echo "$message" | tee -a "$log_file"
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

# Function to check the status of the last executed command
check_status() {
    local status=$1
    if [ $status -ne 0 ]; then
        log "Error: Last command failed with status $status" false
        success=false  # Mark script as failed
    fi
}

# Initialize status variable
success=true

# Example usage
log "---------------Script Timeline -------------------------" false
log "Script started"
log "Adding and committing scripts to git"

# Add and commit files to Git
add_and_commit_to_git "Automated commit before running script"
check_status $?  # Check if the Git commit was successful
log "Files added and committed to Git"  # No timestamp for this message

# Get the current Git commit hash
log "Fetching the current git hash"
commit_hash=$(get_git_commit)
log "Current Git commit hash: $commit_hash"  # No timestamp for this message

# Store the commands in variables
skimmer_command="python skimmer.py $skimmer_args -o $output_file"
message_bot_command="python messegeBOT.py"

# Log and execute the commands
log "$skimmer_command"
eval "$skimmer_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the skimmer command was successful

log "$message_bot_command"
eval "$message_bot_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the message bot command was successful

log "Script completed"
log "===========================================================" false

log "---------------Summary -------------------------"
if [ -n "$description" ]; then
    log "Description: $description" false
fi
log "Timestamp: $timestamp" false
log "Log filename: $log_file" false
log "Skimmer output filename: $output_file" false
log "Git hash: $commit_hash" false

if [ "$success" = true ]; then
    log "Script Status: Successfully run" false
else
    log "Script Status: Failed" false
fi

log "===========================================================" false
