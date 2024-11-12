#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Compute suffix array using provided Karkkainen-Sanders implementation (src/KarkSand.py).

# Implement querying (lecture 01 slide 48). For now do O(n * log m).
# Improve to O(n + log m) if possible using https://www.mimuw.edu.pl/~szczurek/TSG2/04_suffix_arrays.pdf
# Be aware of memory usage here!

# Implement index-assisted approximate matching (lecture 03 slide 55).

# Tweak the parameters to achieve maximum efficiency and accuracy.

from KarkSand import direct_kark_sort
from collections import defaultdict
from dataclasses import dataclass
from typing import List
import sys
import os
import numpy
import matplotlib.pyplot as plt

# DP algorithm adapted from Langmead's notebooks
def _trace(DP, pattern, reference):
    print(f"trace: {len(pattern)}, {len(reference)}")
    ''' Backtrace edit-distance matrix D for strings x and y '''
    index_pattern, index_reference = len(pattern), len(reference)

    while index_pattern > 0:
        # Set initial large values for comparisons
        diagonal = vertical = horizontal = sys.maxsize
        if index_pattern > 0 and index_reference > 0:
            # Diagonal movement: match or substitution
            delta = int(pattern[index_pattern - 1] != reference[index_reference - 1])
            diagonal = DP[index_pattern - 1, index_reference - 1] + delta
        if index_pattern > 0:
            # Vertical movement: insertion
            vertical = DP[index_pattern - 1, index_reference] + 1
        if index_reference > 0:
            # Horizontal movement: deletion
            horizontal = DP[index_pattern, index_reference - 1] + 1

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
    print(f"trace completed.")
    return index_reference

def __kEditDp(pattern, reference, edist):
    print(f"kEditDp - pattern: {len(pattern)} | reference: {len(reference)}")
    DP_edist = numpy.zeros((len(pattern) + 1, len(reference) + 2), dtype=int)
    DP_edist[1 : , 0] = range(1, len(pattern) + 1)
    DP_edist[len(pattern), : ] = [len(pattern) + len(reference) for _ in range(0, len(reference) + 2)]

    DP_start = numpy.zeros((len(pattern) + 1, len(reference) + 2), dtype=int)
    DP_start[0, : ] = range(0, len(reference) + 2)

    DP_delta = numpy.array(
        [
            [
                int(pattern[row] != reference[col])
                for col in range(len(reference))
            ]
            for row in range(len(pattern))
        ],
        dtype=int
    )

    print("kEditDp before first loops.")
    for row in range(1, len(pattern) + 1):
        DP_edist[row, max(0, row - edist - 1)] = row
        DP_edist[row, min(len(reference) + 1, len(reference) + len(pattern) + edist + row)] = \
            len(pattern) + len(reference)
        DP_start[row, min(len(reference) + 1, len(reference) + len(pattern) + edist + row)] = \
            min(len(reference), len(reference) + len(pattern) + edist + row - 1)

        for col in range(
            max(1, row - edist),
            min(len(reference) + 1, len(reference) + len(pattern) + edist + row)
        ):
            diagonal_edist = DP_edist[row - 1, col - 1] + DP_delta[row - 1, col - 1]
            vertical_edist = DP_edist[row - 1, col]
            horizontal_edist = DP_edist[row, col - 1]

            next_edist = min(diagonal_edist, vertical_edist, horizontal_edist)
            DP_edist[row, col] = next_edist

            if next_edist == diagonal_edist:
                DP_start[row, col] = DP_start[row - 1, col - 1]
            elif next_edist == vertical_edist:
                DP_start[row, col] = DP_start[row - 1, col]
            else: # next_edist == horizontal_edist
                DP_start[row, col] = DP_start[row, col - 1]

    print("kEditDp first loops completed.")
    print("kEditDp before second loop.")
    best_edist, best_start, best_end = len(pattern) + len(reference), 0, len(reference)
    for col in range(len(reference) + 1):
        if DP_edist[len(pattern), col] < best_edist:
            best_edist, best_start = DP_edist[len(pattern), col], DP_start[len(pattern), col]

    print("kEditDp second loop completed.")
    return best_edist, best_start, best_end

def _kEditDp(pattern, reference, edist):
    print("Start kEditDP.")
    len_pattern, len_reference = len(pattern), len(reference)

    DP = numpy.zeros((len_pattern + 1, len_reference + 1), dtype=int)
    DP[1:, 0] = numpy.arange(1, len_pattern + 1)

    delta = (numpy.array(list(pattern), dtype='<U1')[ : , None] != numpy.array(list(reference), dtype='<U1')).astype(int)

    for diag in range(2, len_pattern + len_reference + 1):
        row_start = min(len_pattern, diag - 1)
        row_end = max(0, diag - len_reference - 1)

        col_start = max(1, diag - len_pattern)
        col_end = min(len_reference + 1, diag)

        row_indexes = numpy.arange(row_start, row_end, -1)
        col_indexes = numpy.arange(col_start, col_end)

        DP[row_indexes, col_indexes] = numpy.minimum.reduce([
            DP[row_indexes - 1, col_indexes - 1] + delta[row_indexes - 1, col_indexes - 1],
            DP[row_indexes - 1, col_indexes] + 1,
            DP[row_indexes, col_indexes - 1] + 1
        ])

    best_end = numpy.argmin(DP[len_pattern, : len_reference + 1])
    best_start = _trace(DP, pattern, reference[ : best_end])
    best_edist = DP[len_pattern, best_end]

    print("End kEditDp.")
    return best_edist, best_start, best_end


def __kEditDp(p, t, edist):
    ''' Find the alignment of p to a substring of t with the fewest edits.
        Return the edit distance and the coordinates of the substring. '''
    print(f"kEditDp - p: {len(p)} | t: {len(t)}")
    D = numpy.zeros((len(p)+1, len(t)+2), dtype=int)
    # Note: First row gets zeros.  First column initialized as usual.
    D[1:, 0] = range(1, len(p)+1)
    D[len(p), : ] = 10**9
    print("kEditDp before first loops.")
    for i in range(1, len(p)+1):
        D[i, max(0, i - edist - 1)] = i
        D[i, min(len(t) + 1, len(t) - len(p) + edist + i)]
        for j in range(max(1, i - edist), min(len(t) + 1, len(t) - len(p) + edist + i)):
            delt = 1 if p[i-1] != t[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    print("kEditDp first loops completed.")
    # Find minimum edit distance in last row
    mnJ, mn = None, len(p) + len(t)
    for j in range(len(t)+1):
        if D[len(p), j] < mn:
            mnJ, mn = j, D[len(p), j]
    print("kEditDp second loop completed.")
    # Backtrace; note: stops as soon as it gets to first row
    off = _trace(D, p, t[:mnJ])
    # Return edit distance and t coordinates of aligned substring
    print("kEditDp completed.")
    return mn, off, mnJ

def partition(p, num_parts):
    ''' Divide p into num_parts non-overlapping parts, returning each part and its starting position '''
    part_size = len(p) // num_parts
    return [(p[i*part_size : (i+1)*part_size], i*part_size) for i in range(num_parts)]

@dataclass(order=True)
class Shingle:
    pattern_occurence_exact_list: List[int]
    pattern_occurence_range_list: List[range]
    pattern_offset: int
    shingle: str

    def __init__(self, pattern_offset, shingle):
        self.pattern_occurence_exact_list = []
        self.pattern_occurence_range_list = []
        self.pattern_offset = pattern_offset
        self.shingle = shingle

    def __iter__(self):
        yield self.pattern_occurence_exact_list
        yield self.pattern_occurence_range_list
        yield self.pattern_offset
        yield self.shingle

    def __hash__(self):
        return hash((self.pattern_offset, self.shingle))

    def updatePatternOccurenceRange(self, edist, reference_occurence_index_list, reference_length):
        self.pattern_occurence_exact_list = [
            index - self.pattern_offset
            for index in reference_occurence_index_list
        ]
        self.pattern_occurence_range_list = [
            range(
                max(0, index - self.pattern_offset - edist),
                min(reference_length, index - self.pattern_offset + edist + 1)
            )
            for index in reference_occurence_index_list
        ]

class SimpleIndex:
    def __init__(self, reference):
        self.reference = reference
        print("Initializing index...")
        self.suffix_array = direct_kark_sort(reference)
        print("Initialization completed.")

    def __binarySearchSuffixArray(self, pattern):
        assert self.reference[-1] == '$'  # t already has terminator
        assert len(self.reference) == len(self.suffix_array)  # sa is the suffix array for t
        if len(self.reference) == 1:
            return 1

        l, r = 0, len(self.suffix_array)  # invariant: sa[l] < p < sa[r]

        while True:
            c = (l + r) // 2
            # determine whether p < T[sa[c]:] by doing comparisons
            # starting from left-hand sides of p and T[sa[c]:]
            plt = True  # assume p < T[sa[c]:] until proven otherwise
            i = 0

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
            else:
                if c == r - 1:
                    return r
                l = c


    # Returns a list of shingles obtained from pattern, with edit distances of 0 and 1
    def __shingles(self, pattern, edist):
        __SHINGLE_STEP = 10
        __SHINGLE_LENGTH_MIN = 10
        __NUCLEOTIDES = "GACT"

        # Secure the shingling parameters to avoid weird behavior
        shingle_length = max(2 * (len(pattern) // edist), __SHINGLE_LENGTH_MIN)
        shingle_step = min(__SHINGLE_STEP, shingle_length)

        shingle_set = set()
        neighbor_shingle_set = set()

        # Compute shingles with expected edit distance = 0
        for i in range(0, len(pattern) - shingle_length + 1, shingle_step):
            shingle_set.add(Shingle(i, pattern[i : i + shingle_length]))

        # If the last shingle does not reach the end of pattern, add an extra one
        if pattern[-shingle_length : ] not in shingle_set:
            shingle_set.add(Shingle(len(pattern) - shingle_length, pattern[-shingle_length : ]))

        # Compute shingles with expected edit distance = 1
        for _, _, offset, shingle in shingle_set:
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


    # Returns a list of positions at which the pattern occurs in reference
    def __querySuffixArray(self, pattern):
        # Helper for incrementing the ASCII code of last character by 1
        def __increment_last_character(string):
            return string[ : -1] + chr(ord(string[-1]) + 1)

        matching_range = range(
            self.__binarySearchSuffixArray(pattern),
            self.__binarySearchSuffixArray(__increment_last_character(pattern)))

        return [self.suffix_array[i] for i in matching_range]

    def query(self, pattern, edist, read_id):
        shingle_list = self.__shingles(pattern, edist)

        occurence_map = defaultdict(int)

        for shingle in shingle_list:
            occurence_list = self.__querySuffixArray(shingle.shingle)
            shingle.updatePatternOccurenceRange(edist, occurence_list, len(self.reference))
            for occurence_range in shingle.pattern_occurence_range_list:
                occurence_map[occurence_range.start] += 1
                occurence_map[occurence_range.stop] -= 1

        cumulative_occurence_count_list = []
        cumulative_occurence_count = 0
        cumulative_occurence_max_count = 0
        cumulative_occurence_max_index = None

        for index, count in sorted(occurence_map.items()):
            cumulative_occurence_count += count
            cumulative_occurence_count_list.append((index, cumulative_occurence_count))

            if cumulative_occurence_count > cumulative_occurence_max_count:
                cumulative_occurence_max_count = cumulative_occurence_count
                cumulative_occurence_max_index = index

        reference_substring_start = max(0, cumulative_occurence_max_index - 2 * edist)
        reference_substring_end = min(len(self.reference), cumulative_occurence_max_index + len(pattern) + 2 * edist)
        reference_substring = self.reference[reference_substring_start : reference_substring_end]

        # return occurence_map
        # self.__plot_occurence_map(cumulative_occurence_count_list, read_id)

        best_edist, occurence_start, occurence_end = _kEditDp(pattern, reference_substring, 2 * edist)
        return [(
            best_edist,
            reference_substring_start + occurence_start,
            reference_substring_start + occurence_end
        )]

    def __plot_occurence_map(self, cumulative_occurence_count_list, read_id):
        index_list = [index for index, _ in cumulative_occurence_count_list]
        count_list = [count for _, count in cumulative_occurence_count_list]

        plt.figure(figsize=(12, 6))
        plt.bar(index_list, count_list, width=1, color='lightblue', edgecolor='black')
        plt.xlabel("Reference Index")
        plt.ylabel("Occurrence Count")
        plt.title("Shingle Occurrence Counts by Reference Index")

        output_directory = "tests/large/output/plots"
        os.makedirs(output_directory, exist_ok=True)

        plt.tight_layout()
        plt.savefig(f"{output_directory}/plot{read_id}.png")
        plt.close()


from Bio import SeqIO
from sys import argv

ERROR_RATE = 0.1
REFERENCE_SUFFIX = "$"

seq_rec_list = [seq_record for seq_record in SeqIO.parse(argv[1], "fasta")]
index = SimpleIndex(str(seq_rec_list[0].seq) + REFERENCE_SUFFIX)
del seq_rec_list

fout = open(argv[3], "w")
reads = SeqIO.parse(argv[2], "fasta")
for read in reads:
    print(f"Processing read: {read.id}")
    hits = index.query(str(read.seq), int(len(read.seq) * ERROR_RATE), read.id)
    if hits:
        fout.write("{}\t{}\t{}\n".format(read.id, hits[0][1], hits[0][2]))
fout.close()
