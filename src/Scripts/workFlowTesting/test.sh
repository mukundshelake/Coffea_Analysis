#!/bin/bash

# Generate a timestamp for the log filename
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")
log_file="script_log_$timestamp.log"



# Function to log messages with an optional timestamp
log() {
    local message=$1
    local with_timestamp=${2:-true}  # Default to true if not specified

    if [ "$with_timestamp" = true ]; then
        local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
        echo "$timestamp - $message" >> "$log_file"
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

# Your script commands
log "Running command: ls -l"
ls -l | tee -a "$log_file"  # Capture and log output

log "Script completed"
