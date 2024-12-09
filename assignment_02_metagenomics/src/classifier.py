import os
import sys
import gzip
import time
from collections import defaultdict
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


class Classifier:
    def __init__(
            self,
            training_file,
            testing_file,
            output_file,
            k=7,
            downsample_percentile=65):
        self.training_file = training_file
        self.testing_file = testing_file
        self.output_file = output_file
        self.fasta_dir = os.path.dirname(training_file)  # Automatically infer path
        self.k = k
        self.downsample_percentile = downsample_percentile
        self.training_data = None
        self.testing_data = None
        self.classes = []
        self.class_kmer_profiles = {}
        self.total_reads = 0  # Track the total number of reads processed

    def load_data(self):
        """Load training and testing data from TSV files, considering only the first two columns."""
        try:
            # Parse training data with only the first two columns
            print("Loading training data...")
            self.training_data = pd.read_csv(self.training_file, sep="\t", usecols=[0, 1], names=["fasta_file", "geo_loc_name"], skiprows=1)
            print("Training Data Loaded Successfully!")

            # Parse testing data with only the fasta_file column
            print("Loading testing data...")
            self.testing_data = pd.read_csv(self.testing_file, sep="\t", usecols=[0], names=["fasta_file"], skiprows=1)
            print("Testing Data Loaded Successfully!")
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)

    def extract_classes(self):
        """Extract unique classes (geo_loc_name) from training data."""
        self.classes = self.training_data["geo_loc_name"].unique()
        print(f"Classes Extracted: {self.classes}")


    def parse_fasta(self, fasta_file):
        """Read and extract sequences from a gzipped FASTA file."""
        sequences = []
        with gzip.open(fasta_file, "rt") as f:
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                sequences.append(sequence)
        self.total_reads += len(sequences)  # Update total reads count
        return sequences

    def generate_kmers(self, sequence, k):
        """Generate k-mers from a single sequence."""
        return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

    def compute_kmer_profile(self, sequences):
        """
        Compute k-mer counts for a list of sequences, with optional downsampling based on a percentile.
        """
        kmer_counts = defaultdict(int)

        # Count k-mers
        for seq in sequences:
            for kmer in self.generate_kmers(seq, self.k):
                kmer_counts[kmer] += 1

        # Downsample very frequent k-mers based on self.downsample_percentile
        # Extract all k-mer counts and calculate the downsample threshold
        counts = np.array(list(kmer_counts.values()))
        threshold = np.percentile(counts, self.downsample_percentile)

        # Apply the downsampling threshold
        for kmer in kmer_counts:
            if kmer_counts[kmer] > threshold:
                kmer_counts[kmer] = threshold

        return kmer_counts

    def aggregate_class_profiles(self):
        """Aggregate k-mer profiles for each class."""
        for class_name in self.classes:
            class_files = self.training_data[self.training_data['geo_loc_name'] == class_name]['fasta_file']
            aggregated_kmer_counts = defaultdict(int)
            for fasta_file in class_files:
                full_path = os.path.join(self.fasta_dir, fasta_file)
                sequences = self.parse_fasta(full_path)
                class_kmers = self.compute_kmer_profile(sequences)
                for kmer, count in class_kmers.items():
                    aggregated_kmer_counts[kmer] += count
            self.class_kmer_profiles[class_name] = aggregated_kmer_counts
        print("Aggregated k-mer profiles for all classes.")

    def generate_likelihoods(self):
        """Generate likelihoods for each test sample based on k-mer similarity."""
        results = []
        for fasta_file in self.testing_data['fasta_file']:
            full_path = os.path.join(self.fasta_dir, fasta_file)
            sequences = self.parse_fasta(full_path)
            test_kmer_counts = self.compute_kmer_profile(sequences)

            # Create a universal k-mer set
            all_kmers = set(test_kmer_counts.keys())
            for class_kmers in self.class_kmer_profiles.values():
                all_kmers.update(class_kmers.keys())

            # Build aligned vectors
            test_vector = np.array([test_kmer_counts.get(kmer, 0) for kmer in all_kmers])
            class_scores = {}

            for class_name, class_kmers in self.class_kmer_profiles.items():
                class_vector = np.array([class_kmers.get(kmer, 0) for kmer in all_kmers])

                # Calculate similarity (cosine similarity)
                class_scores[class_name] = cosine_similarity([test_vector], [class_vector])[0][0]

            # Normalize scores
            total_score = sum(class_scores.values())
            normalized_scores = [class_scores[class_name] / total_score for class_name in self.classes]

            results.append([fasta_file] + normalized_scores)

        return results


    def write_output(self, results):
        """Write the classifier output to a TSV file."""
        output_columns = ['fasta_file'] + list(self.classes)
        output_df = pd.DataFrame(results, columns=output_columns)
        output_df.to_csv(self.output_file, sep="\t", index=False)
        print(f"Output written to {self.output_file}")

    def run(self):
        """Run the classification pipeline and measure performance."""
        start_time = time.time()

        self.load_data()
        self.extract_classes()
        self.aggregate_class_profiles()
        results = self.generate_likelihoods()
        self.write_output(results)

        end_time = time.time()
        elapsed_time = end_time - start_time  # Total runtime in seconds

                # Calculate elapsed time in minutes and seconds
        minutes, seconds = divmod(elapsed_time, 60)
        formatted_elapsed_time = f"{int(minutes)}m {int(seconds)}s"

        # Calculate runtime per 1M reads
        runtime_per_million_reads = elapsed_time / (self.total_reads / 1_000_000)
        minutes, seconds = divmod(runtime_per_million_reads, 60)
        formatted_runtime_per_million_reads = f"{int(minutes)}m {int(seconds)}s"

        # Print the results
        print(f"Total reads processed: {self.total_reads}")
        print(f"Runtime: {formatted_elapsed_time}")
        print(f"Runtime per 1M reads: {formatted_runtime_per_million_reads}")


def main():
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python3 classifier.py training_data.tsv testing_data.tsv output.tsv")
        sys.exit(1)

    training_file = sys.argv[1]
    testing_file = sys.argv[2]
    output_file = sys.argv[3]

    # Run the classifier
    classifier = Classifier(training_file, testing_file, output_file)
    classifier.run()


if __name__ == "__main__":
    main()
