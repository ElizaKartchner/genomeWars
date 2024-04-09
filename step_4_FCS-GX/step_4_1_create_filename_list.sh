#!/bin/bash

# Check if an argument is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <output_file>"
    echo "  - directory: The directory containing files"
    echo "  - output_file: Where to write the filenames"
    exit 1
fi

# Get the directory path from the command line argument
dir="$1"
out="$2"

# List all files in the directory and save to a file
echo "Creating list of filenames"
ls -1 $dir > $out
echo "Filename list created"
