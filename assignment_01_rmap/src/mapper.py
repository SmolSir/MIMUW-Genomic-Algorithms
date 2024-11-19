#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Compute suffix array using provided Karkkainen-Sanders implementation (src/KarkSand.py).

# Implement querying (lecture 01 slide 48). For now do O(n * log m).
# Improve to O(n + log m) if possible using https://www.mimuw.edu.pl/~szczurek/TSG2/04_suffix_arrays.pdf
# Be aware of memory usage here!

# Implement index-assisted approximate matching (lecture 03 slide 55).

# Tweak the parameters to achieve maximum efficiency and accuracy.

from Bio import SeqIO
from collections import defaultdict
from dataclasses import dataclass
from KarkSand import direct_kark_sort
from sys import argv, maxsize
from time import time
from typing import List

import math
import numpy
import os


class KEditDP:
    ''' Class for finding position of the closest edit-distance-wise occurence
    of pattern within reference '''
    def __init__(self, pattern, reference, edist):
        self.pattern = pattern
        self.reference = reference
        self.edist = edist
        self.__initMatrix()

    def __initMatrix(self):
        ''' Initialise the DP matrix's first row and column, as well as the delta matrix '''
        len_pattern, len_reference = len(self.pattern), len(self.reference)

        self.DP = numpy.zeros((len_pattern + 1, len_reference + 1), dtype=int)
        self.DP[1:, 0] = numpy.arange(1, len_pattern + 1)

        self.DP_delta = (
            numpy.array(list(self.pattern), dtype='<U1')[ : , None] !=
            numpy.array(list(self.reference), dtype='<U1')
        ).astype(bool)

    def __trace(self, best_reference_occurence_index):
        ''' Backtrace edit-distance matrix D for strings x and y '''

        reference_prefix = self.reference[ : best_reference_occurence_index]
        index_pattern, index_reference = len(self.pattern), len(reference_prefix)

        while index_pattern > 0:
            # Set initial large values for comparisons
            diagonal = vertical = horizontal = maxsize
            if index_pattern > 0 and index_reference > 0:
                # Diagonal movement: match or substitution
                delta = int(self.pattern[index_pattern - 1] != reference_prefix[index_reference - 1])
                diagonal = self.DP[index_pattern - 1, index_reference - 1] + delta
            if index_pattern > 0:
                # Vertical movement: insertion
                vertical = self.DP[index_pattern - 1, index_reference] + 1
            if index_reference > 0:
                # Horizontal movement: deletion
                horizontal = self.DP[index_pattern, index_reference - 1] + 1

            # Determine the direction to move
            if diagonal <= vertical and diagonal <= horizontal:
                # Diagonal was best
                index_pattern -= 1; index_reference -= 1
            elif vertical <= horizontal:
                # Vertical was best
                index_pattern -= 1
            else:
                # Horizontal was best
                index_reference -= 1

        # Return the index offset for alignment start in the reference
        return index_reference

    def compute(self):
        ''' Computes the entire DP matrix and traces down the best edit-distance-wise occurence '''
        len_pattern, len_reference = len(self.pattern), len(self.reference)

        # Compute in diagonals to take advantage of numpy vectorisations
        for diag in range(2, len_pattern + len_reference + 1):
            row_start = min(len_pattern, diag - 1)
            row_end = max(0, diag - len_reference - 1)

            col_start = max(1, diag - len_pattern)
            col_end = min(len_reference + 1, diag)

            row_indexes = numpy.arange(row_start, row_end, -1)
            col_indexes = numpy.arange(col_start, col_end)

            self.DP[row_indexes, col_indexes] = numpy.minimum.reduce([
                self.DP[row_indexes - 1, col_indexes - 1] + self.DP_delta[row_indexes - 1, col_indexes - 1],
                self.DP[row_indexes - 1, col_indexes] + 1,
                self.DP[row_indexes, col_indexes - 1] + 1
            ])

        best_end = numpy.argmin(self.DP[len_pattern, : len_reference + 1])
        best_start = self.__trace(best_end)
        best_edist = self.DP[len_pattern, best_end]

        return best_edist, best_start, best_end


@dataclass(order=True)
class Shingle:
    ''' Class representing substrings of pattern that is being matched against the reference '''
    pattern_occurence_range_list: List[range]
    pattern_offset: int
    shingle: str

    def __init__(self, pattern_offset, shingle):
        self.pattern_occurence_range_list = []
        self.pattern_offset = pattern_offset
        self.shingle = shingle

    def __iter__(self):
        yield self.pattern_occurence_range_list
        yield self.pattern_offset
        yield self.shingle

    def __hash__(self):
        return hash((self.pattern_offset, self.shingle))

    def updatePatternOccurenceRange(self, edist, reference_occurence_index_list, reference_length):
        ''' Given a list of shingle occurence positions in reference and maximal expected
        edit distance, for each occurence computes a range representing the possible
        starting positions of the shingle's pattern in reference '''
        self.pattern_occurence_range_list = [
            range(
                max(0, index - self.pattern_offset - edist),
                min(reference_length, index - self.pattern_offset + edist + 1)
            )
            for index in reference_occurence_index_list
        ]


class Index:
    ''' Class responsible for computing queries '''
    def __init__(self, reference):
        print("Initialising index...")
        index_start_time = time()

        self.reference = reference
        self.suffix_array = direct_kark_sort(reference)

        index_end_time = time()
        index_time = index_end_time - index_start_time
        print(f"Initialisation completed in: {int(index_time // 60)}m{(index_time % 60):.3f}s")

    def __binarySearchSuffixArray(self, pattern):
        ''' Computes the leftmost index where pattern could be inserted while keeping
        sorted order '''
        assert self.reference[-1] == '$'  # t already has terminator
        assert len(self.reference) == len(self.suffix_array)  # sa is the suffix array for t
        if len(self.reference) == 1:
            return 1

        l, r = 0, len(self.suffix_array)  # invariant: sa[l] < p < sa[r]
        lcp_l, lcp_r = 0, 0

        while True:
            c = (l + r) // 2
            # determine whether p < T[sa[c]:] by doing comparisons
            # starting from left-hand sides of p and T[sa[c]:]
            plt = True  # assume p < T[sa[c]:] until proven otherwise
            i = min(lcp_l, lcp_r)

            while i < len(pattern) and self.suffix_array[c] + i < len(self.reference):
                if pattern[i] < self.reference[self.suffix_array[c] + i]:
                    break  # p < T[sa[c]:]
                elif pattern[i] > self.reference[self.suffix_array[c] + i]:
                    plt = False
                    break  # p > T[sa[c]:]
                i += 1  # tied so far

            if plt:
                if c == l + 1:
                    return c
                r = c
                lcp_r = i
            else:
                if c == r - 1:
                    return r
                l = c
                lcp_l = i

    def __shingles(self, pattern, edist):
        ''' Computes a list of shingles obtained from pattern, with edit distances of 0 and 1 '''
        __SHINGLE_STEP = 10
        __SHINGLE_LENGTH_MIN = 10
        __NUCLEOTIDES = "GACT"

        # Secure the shingling parameters to avoid weird behavior
        shingle_length = max(2 * (len(pattern) // (edist + 1)), __SHINGLE_LENGTH_MIN)
        shingle_step = min(__SHINGLE_STEP, shingle_length // 2)

        shingle_set = set()
        neighbor_shingle_set = set()

        # Compute shingles with expected edit distance = 0
        for i in range(0, len(pattern) - shingle_length + 1, shingle_step):
            shingle_set.add(Shingle(i, pattern[i : i + shingle_length]))

        # If the last shingle does not reach the end of pattern, add an extra one
        if pattern[-shingle_length : ] not in shingle_set:
            shingle_set.add(Shingle(len(pattern) - shingle_length, pattern[-shingle_length : ]))

        # Compute shingles with expected edit distance = 1
        for _, offset, shingle in shingle_set:
            # Compute shingles with one substitution
            for i in range(len(shingle)):
                for nucleotide in __NUCLEOTIDES:
                    if shingle[i] != nucleotide:
                        neighbor_shingle = shingle[ : i] + nucleotide + shingle[i + 1 : ]
                        neighbor_shingle_set.add(Shingle(offset, neighbor_shingle))

            # Compute shingles with one insertion
            for i in range(len(shingle) + 1):
                for nucleotide in __NUCLEOTIDES:
                    neighbor_shingle = shingle[ : i] + nucleotide + shingle[i : ]
                    neighbor_shingle_set.add(Shingle(offset, neighbor_shingle))

            # Compute shingles with one deletion
            for i in range(len(shingle)):
                neighbor_shingle = shingle[ : i] + shingle[i + 1 : ]
                neighbor_shingle_set.add(Shingle(offset, neighbor_shingle))

        # Merge both sets and return as list
        shingle_set.update(neighbor_shingle_set)
        return list(shingle_set)

    def __querySuffixArray(self, pattern):
        ''' Computes a list of pattern's occurence positions in reference '''
        increment_pattern = pattern[ : -1] + chr(ord(pattern[-1]) + 1)
        matching_range = range(
            self.__binarySearchSuffixArray(pattern),
            self.__binarySearchSuffixArray(increment_pattern))

        return [self.suffix_array[i] for i in matching_range]

    def query(self, pattern, edist, read_id):
        ''' Finds the closest matching occurence edit-distance-wise of pattern in reference '''
        print(f"Processing read: {read_id}")

        shingle_list = self.__shingles(pattern, edist)
        occurence_map = defaultdict(int)

        for shingle in shingle_list:
            occurence_list = self.__querySuffixArray(shingle.shingle)
            shingle.updatePatternOccurenceRange(edist, occurence_list, len(self.reference))
            for occurence_range in shingle.pattern_occurence_range_list:
                occurence_map[occurence_range.start] += 1
                occurence_map[occurence_range.stop] -= 1

        cumulative_occurence_count = 0
        cumulative_occurence_max_count = 0
        cumulative_occurence_max_index = None

        for index, count in sorted(occurence_map.items()):
            cumulative_occurence_count += count
            if cumulative_occurence_count > cumulative_occurence_max_count:
                cumulative_occurence_max_count = cumulative_occurence_count
                cumulative_occurence_max_index = index

        reference_substring_start = max(0, cumulative_occurence_max_index - 2 * edist)
        reference_substring_end = min(len(self.reference), cumulative_occurence_max_index + len(pattern) + 2 * edist)
        reference_substring = self.reference[reference_substring_start : reference_substring_end]

        # self.__plot_occurence_map(occurence_map, read_id)

        kEditDP = KEditDP(pattern, reference_substring, 2 * edist)
        best_edist, occurence_start, occurence_end = kEditDP.compute()

        return Result(
            read_id,
            best_edist,
            reference_substring_start + occurence_start,
            reference_substring_start + occurence_end
        )

    def __plot_occurence_map(self, occurence_map, read_id):
        ''' Plots the expected occurence start positions in reference '''
        import matplotlib.pyplot as plt

        cumulative_occurence_count_list = []
        cumulative_occurence_count = 0

        for index, count in sorted(occurence_map.items()):
            cumulative_occurence_count += count
            cumulative_occurence_count_list.append((index, cumulative_occurence_count))

        index_list = [index for index, _ in cumulative_occurence_count_list]
        count_list = [count for _, count in cumulative_occurence_count_list]

        plt.figure(figsize=(12, 6))
        plt.bar(index_list, count_list, width=1, color='lightblue', edgecolor='black')
        plt.xlabel("Reference Index")
        plt.ylabel("Occurrence Count")
        plt.title("Shingle Occurrence Counts by Reference Index")

        output_directory = "tests/plots"
        os.makedirs(output_directory, exist_ok=True)

        plt.tight_layout()
        plt.savefig(f"{output_directory}/plot{read_id}.png")
        plt.close()


class Result:
    ''' Class representing the results of a single query '''
    def __init__(self, read_id, best_edist, best_occurence_start, best_occurence_end):
        self.read_id = read_id
        self.best_edist = best_edist
        self.best_occurence_start = best_occurence_start
        self.best_occurence_end = best_occurence_end

    def __iter__(self):
        yield self.read_id
        yield self.best_edist
        yield self.best_occurence_start
        yield self.best_occurence_end


ERROR_RATE = 0.1
READ_RELIABLE_THRESHOLD = 0.99
REFERENCE_SUFFIX = "$"

seq_rec_list = [seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
index = Index(str(seq_rec_list[0].seq) + REFERENCE_SUFFIX)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")

results = [index.query(str(read.seq), int(len(read.seq) * ERROR_RATE), read.id) for read in reads]
results.sort(key=lambda result: result.best_edist)

reliable_results = results[ : math.ceil(len(results) * READ_RELIABLE_THRESHOLD)]
reliable_results.sort(key=lambda result: result.read_id)

for read_id, _, best_occurence_start, best_occurence_end in reliable_results:
    fout.write("{}\t{}\t{}\n".format(read_id, best_occurence_start, best_occurence_end))

fout.close()
