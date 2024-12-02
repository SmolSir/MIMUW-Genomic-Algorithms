import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from classifier import Classifier


def main():
    # Parse command-line arguments
    if len(sys.argv) != 5:
        print("Usage: python3 experiment_percentile.py training_data.tsv testing_data.tsv ground_truth.tsv output_dir")
        sys.exit(1)

    training_file = sys.argv[1]
    testing_file = sys.argv[2]
    ground_truth_file = sys.argv[3]
    output_dir = sys.argv[4]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define range for percentile values
    percentile_values = np.arange(1, 101, 1)

    auc_results = {}

    for percentile in percentile_values:
        print(f"Running experiment for percentile={percentile}")
        output_file = os.path.join(output_dir, f"output_p{percentile}.tsv")

        # Run classifier with specified percentile
        classifier = Classifier(training_file, testing_file, output_file)
        classifier.downsample_percentile = percentile
        classifier.run()

        # Evaluate AUC-ROC using evaluate.py
        auc_roc_scores = evaluate_auc(output_file, ground_truth_file)
        auc_results[percentile] = auc_roc_scores
        print(f"AUC-ROC for percentile={percentile}: {auc_roc_scores}")

    # Plot results
    plot_auc_results(auc_results, output_dir)


def evaluate_auc(output_file, ground_truth_file):
    """
    Call evaluate.py script to calculate AUC-ROC scores.
    """
    try:
        command = ["python3", "src/evaluate.py", output_file, ground_truth_file]
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        scores = {}
        for line in result.stdout.splitlines():
            if line.startswith("AUC-ROC for class"):
                parts = line.split(":")
                class_name = parts[0].split()[-1]
                score = float(parts[1].strip())
                scores[class_name] = score
            elif line.startswith("Average AUC-ROC across all classes"):
                avg_score = float(line.split(":")[1].strip())
                scores["Average"] = avg_score

        return scores
    except subprocess.CalledProcessError as e:
        print(f"Error running evaluate.py: {e.stderr}")
        sys.exit(1)


def plot_auc_results(auc_results, output_dir):
    """
    Plot AUC-ROC results for different percentiles.
    """
    percentiles = sorted(auc_results.keys())
    average_scores = [auc_results[p]["Average"] for p in percentiles]

    plt.figure(figsize=(10, 6))
    plt.plot(percentiles, average_scores, marker='o', label="Average AUC-ROC")
    plt.title("AUC-ROC vs. Percentile")
    plt.xlabel("Percentile")
    plt.ylabel("Average AUC-ROC")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "auc_roc_vs_percentile.png"))
    print(f"Plot saved to {os.path.join(output_dir, 'auc_roc_vs_percentile.png')}")


if __name__ == "__main__":
    main()
