#!/bin/bash

# Check if a directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <test_directory>"
    exit 1
fi

# Set up paths based on the provided test directory
TEST_DIR="$1"
OUTPUT_DIR="${TEST_DIR%/}/output"
EXPECTED_DIR="${TEST_DIR%/}/expected"
METRICS_FILE="${TEST_DIR%/}/metrics.txt"

# Check if output and expected directories exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory $OUTPUT_DIR does not exist"
    exit 1
fi

if [ ! -d "$EXPECTED_DIR" ]; then
    echo "Error: Expected directory $EXPECTED_DIR does not exist"
    exit 1
fi

# Memory limit (1GB) and initialize counters
ulimit -v $((1024 * 1024))  # 1GB in KB
total_expected_reads=0
correctly_mapped=0
total_inaccuracies=0
total_output_reads=0

# Function to compare files, allowing ≤20bp coordinate difference
compare_files() {
    local output_file=$1
    local expected_file=$2

    # Read expected results into an associative array
    declare -A expected_set
    while read -r line; do
        read_id=$(echo "$line" | awk '{print $1}')
        start=$(echo "$line" | awk '{print $2}')
        end=$(echo "$line" | awk '{print $3}')
        expected_set["$read_id"]="$start $end"
        ((total_expected_reads++))
    done < "$expected_file"

    file_inaccuracies=0

    # Compare each line in the output file
    while read -r line; do
        ((total_output_reads++))
        read_id=$(echo "$line" | awk '{print $1}')
        output_start=$(echo "$line" | awk '{print $2}')
        output_end=$(echo "$line" | awk '{print $3}')

        # Check if this read exists in the expected output
        if [[ -n "${expected_set[$read_id]}" ]]; then
            expected_start=$(echo "${expected_set[$read_id]}" | awk '{print $1}')
            expected_end=$(echo "${expected_set[$read_id]}" | awk '{print $2}')

            # Calculate absolute differences
            start_diff=$(( output_start > expected_start ? output_start - expected_start : expected_start - output_start ))
            end_diff=$(( output_end > expected_end ? output_end - expected_end : expected_end - output_end ))

            # Check if both differences are within the 20bp threshold
            if (( start_diff <= 20 && end_diff <= 20 )); then
                ((correctly_mapped++))
            else
                ((file_inaccuracies++))
            fi
        else
            ((file_inaccuracies++))
        fi
    done < "$output_file"

    total_inaccuracies=$((total_inaccuracies + file_inaccuracies))
}

# Loop through each output file and compare with expected
for output_file in "$OUTPUT_DIR"/*.out; do
    base_name=$(basename "$output_file" .out)
    expected_file="$EXPECTED_DIR/$base_name.txt"

    if [ -f "$expected_file" ]; then
        echo "Comparing $output_file with $expected_file"
        compare_files "$output_file" "$expected_file"
    else
        echo "Warning: Expected file $expected_file not found for $output_file"
    fi
done

# Calculate final metrics
mapped_percentage=$(echo "scale=2; ($correctly_mapped / $total_expected_reads) * 100" | bc)
inaccuracy_percentage=$(echo "scale=2; ($total_inaccuracies / $total_output_reads) * 100" | bc)

# Output final metrics to both console and metrics.txt
{
    echo "Final Metrics:"
    echo "Total Reads Mapped: $correctly_mapped / $total_expected_reads ($mapped_percentage%)"
    echo "Total Inaccuracies: $total_inaccuracies ($inaccuracy_percentage%)"

    # Check assignment requirements
    if (( $(echo "$mapped_percentage >= 80" | bc -l) )); then
        echo "Mapped Percentage meets the requirement (≥80%)"
    else
        echo "Mapped Percentage does not meet the requirement (≥80%)"
    fi

    if (( $(echo "$inaccuracy_percentage <= 1" | bc -l) )); then
        echo "Inaccuracy Percentage meets the requirement (≤1%)"
    else
        echo "Inaccuracy Percentage does not meet the requirement (≤1%)"
    fi
} | tee "$METRICS_FILE"  # Output to both console and file

echo "Metrics written to $METRICS_FILE"
