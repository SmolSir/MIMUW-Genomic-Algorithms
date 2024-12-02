import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from classifier import Classifier


def main():
    # Parse command-line arguments
    if len(sys.argv) != 5:
        print("Usage: python3 experiment_filter.py training_data.tsv testing_data.tsv ground_truth.tsv output_dir")
        sys.exit(1)

    training_file = sys.argv[1]
    testing_file = sys.argv[2]
    ground_truth_file = sys.argv[3]
    output_dir = sys.argv[4]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define ranges for threshold and proportion
    threshold_values = np.arange(0.05, 1.0, 0.05)
    proportion_values = np.arange(0.02, 1.0, 0.02)

    auc_results = {}

    for threshold in threshold_values:
        threshold = round(threshold, 2)
        for proportion in proportion_values:
            proportion = round(proportion, 2)
            print(f"Running experiment for threshold={threshold}, proportion={proportion}")
            output_file = os.path.join(output_dir, f"output_t{threshold}_p{proportion}.tsv")

            # Run classifier with threshold and proportion parameters
            classifier = Classifier(training_file, testing_file, output_file)
            classifier.filtering_threshold = threshold
            classifier.filtering_proportion = proportion
            classifier.run()

            # Evaluate AUC-ROC using evaluate.py
            auc_roc_scores = evaluate_auc(output_file, ground_truth_file)
            auc_results[(threshold, proportion)] = auc_roc_scores
            print(f"AUC-ROC for threshold={threshold}, proportion={proportion}: {auc_roc_scores}")

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
    Plot AUC-ROC results for different thresholds and proportions with better visualization.
    """
    thresholds = sorted(set(k[0] for k in auc_results))
    proportions = sorted(set(k[1] for k in auc_results))
    avg_aucs = np.zeros((len(thresholds), len(proportions)))

    for i, threshold in enumerate(thresholds):
        for j, proportion in enumerate(proportions):
            avg_aucs[i, j] = auc_results[(threshold, proportion)]["Average"]

    plt.figure(figsize=(100, 40))
    plt.imshow(avg_aucs, cmap="coolwarm", interpolation="nearest", aspect="auto")
    plt.colorbar(label="Average AUC-ROC")
    plt.xticks(range(len(proportions)), [f"{p:.2f}" for p in proportions], rotation=45)
    plt.yticks(range(len(thresholds)), [f"{t:.2f}" for t in thresholds])
    plt.xlabel("Proportion of k-mers Filtered")
    plt.ylabel("Threshold for Common k-mers")
    plt.title("AUC-ROC Heatmap with Annotations")

    # Annotate the heatmap with values
    for i in range(len(thresholds)):
        for j in range(len(proportions)):
            plt.text(j, i, f"{avg_aucs[i, j]:.3f}", ha="center", va="center", color="black", fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "auc_roc_heatmap_annotated.png"))
    print(f"Heatmap saved to {os.path.join(output_dir, 'auc_roc_heatmap_annotated.png')}")



if __name__ == "__main__":
    main()
