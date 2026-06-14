#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <data_directory>"
    echo "Example: $0 tests/sample"
    exit 1
fi

# Assign parameters
DATA_DIR="$1"
OUTPUT_DIR="$DATA_DIR/experiment_k_values"

# Ensure the data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory '$DATA_DIR' does not exist."
    exit 1
fi

# Clear the output directory if it already exists
if [ -d "$OUTPUT_DIR" ]; then
    echo "Clearing existing output directory: $OUTPUT_DIR"
    rm -rf "$OUTPUT_DIR"
fi

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Define input file paths based on the directory structure
TRAINING_DATA="$DATA_DIR/train_data.tsv"
TESTING_DATA="$DATA_DIR/test_data.tsv"
GROUND_TRUTH="$DATA_DIR/test_ground_truth.tsv"

# Check if all required files exist
if [ ! -f "$TRAINING_DATA" ] || [ ! -f "$TESTING_DATA" ] || [ ! -f "$GROUND_TRUTH" ]; then
    echo "Error: One or more required files are missing in '$DATA_DIR'."
    echo "Expected files: train_data.tsv, test_data.tsv, test_ground_truth.tsv"
    exit 1
fi

# Run the experiment script
python3 src/experiment_k_values.py "$TRAINING_DATA" "$TESTING_DATA" "$GROUND_TRUTH" "$OUTPUT_DIR"
