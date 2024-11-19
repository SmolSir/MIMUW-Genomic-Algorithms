#!/bin/bash

# Set locale to avoid issues with comma as decimal separator
export LC_NUMERIC="C"

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

# Initialize counters
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
mapped_percentage=$(echo "scale=6; ($correctly_mapped / $total_expected_reads) * 100" | bc)
inaccuracy_percentage=$(echo "scale=6; ($total_inaccuracies / $total_output_reads) * 100" | bc)

# Format the calculated percentages to two decimal places for output
mapped_percentage=$(printf "%.2f" "$mapped_percentage")
inaccuracy_percentage=$(printf "%.2f" "$inaccuracy_percentage")

# Global variables for points
correctness_points=0
completeness_points=0
time_points=0

# Calculate points for mapping correctness
if (( $(echo "$inaccuracy_percentage == 0" | bc -l) )); then
    correctness_points=3
elif (( $(echo "$inaccuracy_percentage <= 0.2" | bc -l) )); then
    correctness_points=2
elif (( $(echo "$inaccuracy_percentage <= 0.5" | bc -l) )); then
    correctness_points=1
fi

# Calculate points for mapping completeness
if (( $(echo "$mapped_percentage >= 99" | bc -l) )); then
    completeness_points=3
elif (( $(echo "$mapped_percentage >= 95" | bc -l) )); then
    completeness_points=2
elif (( $(echo "$mapped_percentage >= 90" | bc -l) )); then
    completeness_points=1
fi

# Simulate mapping time points calculation (assuming time and r reads are available)
r=$total_output_reads  # Number of reads
processing_time=2  # Placeholder for example; dynamically measure this if possible
max_allowed_time=$(echo "2 + $r / 10" | bc)  # Calculate allowed time as 2 + r/10

if (( $(echo "$processing_time <= 1" | bc -l) )); then
    time_points=3
elif (( $(echo "$processing_time <= 2" | bc -l) )); then
    time_points=2
elif (( $(echo "$processing_time <= 3" | bc -l) )); then
    time_points=1
fi

# Add 3 additional points to the total
total_points=$((2 + correctness_points + completeness_points + time_points + 3))

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

    if (( $(echo "$processing_time <= $max_allowed_time * 60" | bc -l) )); then
        echo "Processing time meets the requirement ($max_allowed_time minutes)"
    else
        echo "Processing time exceeds max allowed time of $max_allowed_time minutes"
    fi

    echo ""
    echo "Points Breakdown:"
    echo "Mapping Correctness Points: $correctness_points"
    echo "Mapping Completeness Points: $completeness_points"
    echo "Mapping Time Points: $time_points"
    echo "Additional Points: 2"
    echo "Total Points: $total_points / 15 pts"
} | tee "$METRICS_FILE"

echo "Metrics written to $METRICS_FILE"
