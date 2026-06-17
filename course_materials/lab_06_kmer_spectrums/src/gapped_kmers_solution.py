import itertools
from math import comb, sqrt

# we consider gapped kmers characterized by two parameters:
#   l -- whole word length including gaps
#   k -- the number of informative (non-gapped) positions in each word
l = 3
k = 2

# the set of sequences
S = ["AAACCC", "AAAAA", "ACC"]


#
#  Naive calculation of the kernel function (similarity score) explicitly using gapped k-mers
#

def string_to_kmer_list(s, k):
    return [s[i:i + k] for i in range(len(s) - k  + 1)]

def mask_kmer(kmer, positions):
    kmer = list(kmer)
    for pos in positions:
        kmer[pos] = "."
    return "".join(kmer)

def kmer_to_gapped_kmer_list(lmer, k):
    l = len(lmer)
    return [mask_kmer(lmer, positions) for positions in itertools.combinations(range(l), l - k)]

def string_to_gapped_kmer_feature_dict(s, l, k):
    d = {}
    for lmer in string_to_kmer_list(s, l):
        for gapped_lmer in kmer_to_gapped_kmer_list(lmer, k):
            d[gapped_lmer] = d.get(gapped_lmer, 0) + 1
    return d

# inner product of vectors corresponding to gapped k-mer feature dicts
def feature_dict_inner_prod(fs1, fs2):
    res = 0
    for k in fs1:
        if k in fs2:
            res += fs1[k] * fs2[k]
    return res

def naive_inner_prod(s1, s2, l, k):
    fs1 = string_to_gapped_kmer_feature_dict(s1, l, k)
    fs2 = string_to_gapped_kmer_feature_dict(s2, l, k)
    return feature_dict_inner_prod(fs1, fs2)

print(naive_inner_prod(S[0], S[1], l, k)) # test, should be 12


def naive_kernel(s1, s2, l, k):
    fs1 = string_to_gapped_kmer_feature_dict(s1, l, k)
    fs2 = string_to_gapped_kmer_feature_dict(s2, l, k)
    return feature_dict_inner_prod(fs1, fs2) / sqrt(feature_dict_inner_prod(fs1, fs1) * feature_dict_inner_prod(fs2, fs2))

def naive_kernel_matrix(S, l, k):
    return [[naive_kernel(s1, s2, l, k) for s2 in S] for s1 in S]

print(naive_kernel_matrix(S, l, k))


#
#  Less naive calculation of the kernel function (similarity score), using only l-mers
#

def strings_inner_prod(s1, s2, l, k):
    res = 0
    for lmer1 in string_to_kmer_list(s1, l):
        for lmer2 in string_to_kmer_list(s2, l):
            mismatches = sum([nucl1 != nucl2 for nucl1, nucl2 in zip(lmer1, lmer2)])
            res += comb(l - mismatches, k)
    return res

print(strings_inner_prod(S[0], S[1], l, k)) # test, should be 12


def strings_kernel_matrix(S, l, k):
    inner_prod_matrix = [[strings_inner_prod(s1, s2, l, k) for s2 in S] for s1 in S]
    return [[inner_prod_matrix[i][j] / sqrt(inner_prod_matrix[i][i] * inner_prod_matrix[j][j])
        for j in range(len(S))]
            for i in range(len(S))]

print(strings_kernel_matrix(S, l, k))


#
#  Calculation of the kernel function (similarity score) using k-mer tree data structure
#

def kmer_tree(S, l):
    # each node is a list:
    #   node[0]: str -- label (not needed! only for convenience)
    #   node[1]: Dict[str, list] -- storing children
    #   node[2]: Dict[int, int] -- for each sequence index, storing the number of times that l-mer appeared in each sequence
    t = ["", {}, {}]

    for i, s in enumerate(S):
        for lmer in string_to_kmer_list(s, l):
            # start from the root
            node = t
            for d, nucleotide in enumerate(lmer):
                if nucleotide not in node[1]:
                    node[1][nucleotide] = [lmer[0:d + 1], {}, {}]
                node = node[1][nucleotide]
            # node is a leaf
            node[2][i] = node[2].get(i, 0) + 1
    return t

def kmer_tree_to_inner_prod_matrix(lt, l, k):
    inner_prod_matrix = [[0 for s2 in S] for s1 in S]

    def dfs_walker(node, depth, other_node_mismatches_list):
        if depth == l:
            for other_node, mismatches in other_node_mismatches_list:
                for i, i_count in node[2].items():
                    for j, j_count in other_node[2].items():
                        if mismatches <= l - k:
                            # print([node[0], i, i_count, other_node[0], j, j_count, mismatches])
                            inner_prod_matrix[i][j] += i_count * j_count * comb(l - mismatches, k)
        else:
            for nucleotide, child_node in node[1].items():
                dfs_walker(child_node, depth + 1,
                    [(other_child_node, mismatches + (0 if other_nucleotide == nucleotide else 1))
                        for other_node, mismatches in other_node_mismatches_list
                            for other_nucleotide, other_child_node in other_node[1].items()])

    dfs_walker(lt, 0, [(lt, 0)])
    return inner_prod_matrix

def kmer_tree_inner_prod_matrix(S, l, k):
    lt = kmer_tree(S, l)
    return kmer_tree_to_inner_prod_matrix(lt, l, k)

print(kmer_tree_inner_prod_matrix(S, l, k)[0][1]) # test, should be 12


def kmer_tree_kernel_matrix(S, l, k):
    lt = kmer_tree(S, l)
    # print(lt)
    inner_prod_matrix = kmer_tree_to_inner_prod_matrix(lt, l, k)
    return [[inner_prod_matrix[i][j] / sqrt(inner_prod_matrix[i][i] * inner_prod_matrix[j][j])
        for j in range(len(S))]
            for i in range(len(S))]

print(kmer_tree_kernel_matrix(S, l, k))
