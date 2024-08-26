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
    shift  # Shift the arguments to the left so that $2 becomes $1, $3 becomes $2, etc.
fi

# Check if a description is provided
if [ -n "$1" ]; then
    description=$1
fi

# Generate a timestamp for the log filename and output filename
timestamp=$(date +"%Y%m%d_%H%M%S${sample_indicator}")
log_file="logs/log_${timestamp}.log"

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
LHSskimmer_command="python LHSskimmer.py $skimmer_args -t $timestamp -c ttbar_SemiLeptonic"
# RHSskimmer_command="python RHSskimmer.py $skimmer_args -t $timestamp -c ttbar_SemiLeptonic"
message_bot_command="python messegeBOT.py"
# LHSextractor_command="python LHSextractor.py -t $timestamp"
# RHSextractor_command="python RHSextractor.py -t $timestamp"
# LHSplotter_command="python LHSplotter.py -t $timestamp"
# RHSplotter_command="python RHSplotter.py -t $timestamp"


# Log and execute the commands
log "$LHSskimmer_command"
eval "$LHSskimmer_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the skimmer command was successful


log "$RHSskimmer_command"
eval "$RHSskimmer_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the skimmer command was successful


log "$message_bot_command"
eval "$message_bot_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the message bot command was successful

log "$LHSextractor_command"
eval "$LHSextractor_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the LHSextractor_command was successful

log "$RHSextractor_command"
eval "$RHSextractor_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the RHSextractor_command was successful

log "$LHSplotter_command"
eval "$LHSplotter_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the LHSplotter_command was successful

log "$RHSplotter_command"
eval "$RHSplotter_command" 2>&1 | tee -a "$log_file"
check_status $?  # Check if the RHSplotter_command was successful


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
