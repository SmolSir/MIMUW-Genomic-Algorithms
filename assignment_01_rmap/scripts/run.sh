#!/bin/bash

# Check if a directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <test_directory>"
    exit 1
fi

# Set up paths based on the provided test directory
TEST_DIR="$1"
SCRIPT_DIR="$(dirname "$0")"
SCRIPT="$SCRIPT_DIR/../src/mapper.py"
INPUT_DIR="${TEST_DIR%/}/input"
OUTPUT_DIR="${TEST_DIR%/}/output"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Check if reference and reads files exist
REFERENCE_FILE=$(find "$TEST_DIR" -name "reference*.fasta" | head -n 1)
if [ -z "$REFERENCE_FILE" ]; then
    echo "Error: No reference file found in $TEST_DIR"
    exit 1
fi

READS_FILES=$(find "$INPUT_DIR" -name "reads*.fasta")
if [ -z "$READS_FILES" ]; then
    echo "Error: No reads files found in $INPUT_DIR"
    exit 1
fi

# Clear the output directory if it exists, or create it if it doesn't
if [ -d "$OUTPUT_DIR" ]; then
    rm -f "$OUTPUT_DIR"/*
else
    mkdir -p "$OUTPUT_DIR"
fi

# Set the memory limit (1GB)
ulimit -v $((1024 * 1024))  # 1GB in KB

# Run tests for each reads file
for READS_FILE in $READS_FILES; do
    # Define output file name based on the reads file name
    OUTPUT_FILE="$OUTPUT_DIR/$(basename "${READS_FILE%.fasta}.out")"

    # Run the mapper.py script with specified arguments
    echo "Running mapper.py on $READS_FILE"
    if ! python3 "$SCRIPT" "$REFERENCE_FILE" "$READS_FILE" "$OUTPUT_FILE"; then
        echo "Error: mapper.py failed on $READS_FILE"
        exit 1
    fi

    echo "Output for $READS_FILE written to $OUTPUT_FILE"
done

echo "All tests in $TEST_DIR completed successfully."
