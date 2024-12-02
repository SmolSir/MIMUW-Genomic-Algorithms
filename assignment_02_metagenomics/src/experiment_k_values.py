import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from classifier import Classifier

def main():
    # Parse command-line arguments
    if len(sys.argv) != 5:
        print("Usage: python3 experiment.py training_data.tsv testing_data.tsv ground_truth.tsv output_dir")
        sys.exit(1)

    training_file = sys.argv[1]
    testing_file = sys.argv[2]
    ground_truth_file = sys.argv[3]
    output_dir = sys.argv[4]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define range of k values to test
    k_values = range(1, 9)  # Experiment with k values from 1 to 10
    auc_results = {}

    for k in k_values:
        print(f"Running experiment for k = {k}")
        output_file = os.path.join(output_dir, f"output_k{k}.tsv")

        # Run classifier
        classifier = Classifier(training_file, testing_file, output_file)
        classifier.k = k
        classifier.run()

        # Use evaluate.py to calculate AUC-ROC
        auc_roc_scores = evaluate_auc(output_file, ground_truth_file)
        auc_results[k] = auc_roc_scores
        print(f"AUC-ROC for k={k}: {auc_roc_scores}")

    # Plot the results
    plot_auc_results(auc_results, output_dir)


def evaluate_auc(output_file, ground_truth_file):
    """
    Call evaluate.py script to calculate AUC-ROC scores.
    """
    try:
        # Construct and execute the evaluate.py command
        command = ["python3", "src/evaluate.py", output_file, ground_truth_file]
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        # Parse AUC scores from the script's output
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
    Plot AUC-ROC results across k values.
    """
    classes = list(next(iter(auc_results.values())).keys())
    classes.remove("Average")  # Exclude "Average" from individual class plots

    plt.figure(figsize=(10, 6))

    for cls in classes:
        class_scores = [auc_results[k][cls] for k in auc_results]
        plt.plot(list(auc_results.keys()), class_scores, marker='o', label=cls)

    # Average AUC-ROC
    average_scores = [auc_results[k]["Average"] for k in auc_results]
    plt.plot(list(auc_results.keys()), average_scores, marker='x', linestyle='--', label="Average", color="black")

    plt.title("AUC-ROC vs. k")
    plt.xlabel("k")
    plt.ylabel("AUC-ROC")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "auc_roc_vs_k.png"))
    print(f"Plot saved to {os.path.join(output_dir, 'auc_roc_vs_k.png')}")


if __name__ == "__main__":
    main()
