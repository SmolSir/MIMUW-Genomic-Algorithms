#!/bin/bash

# Check if a directory parameter is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    echo "Example: $0 tests/sample"
    exit 1
fi

# Assign the directory parameter to a variable
DIR="$1"

# Ensure the directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' does not exist."
    exit 1
fi

# Run the classifier command with the provided directory
python3 src/classifier.py "$DIR/train_data.tsv" "$DIR/test_data.tsv" "$DIR/output/output_data.tsv"
