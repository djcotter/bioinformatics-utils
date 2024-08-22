"""
Deduplicate anchors by clustering them based on their similarity.
Tavor's code
"""

import time
import numpy as np
import scipy
import networkx as nx
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Deduplicate anchors")
    parser.add_argument(
        "--input", type=str, required=True, help="Input file containing anchors"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output file containing deduplicated anchors",
    )
    parser.add_argument(
        "--maxShiftDist",
        type=int,
        default=5,
        help="Maximum number of bases that sequences can be shifted and considered in the same cluster",
    )
    return parser.parse_args()


def bpToInt(x):
    if x == "A":
        return 0
    if x == "T":
        return 1
    if x == "C":
        return 2
    if x == "G":
        return 3
    return 4  # for nan


# i,j-th entry is 1 if i-th row of A equals j-th row of B
def compare_rows_hash(A, B):
    n, k = A.shape
    hash_A = np.array([hash(tuple(row)) for row in A])
    hash_B = np.array([hash(tuple(row)) for row in B])
    return scipy.sparse.csr_matrix(hash_A[:, None] == hash_B)


# main function
# inputs: anchLst, a list of sequences; maxShiftDist, the maximum number of bases that sequences can be shifted and considered in the same cluster
# output: a list of lists, where each list is a cluster of sequences
def clusterAnchors(anchLst, maxShiftDist=5):
    start = time.time()
    bpArr = np.array([[bpToInt(x) for x in s] for s in anchLst], dtype=np.uint8)

    n, k = bpArr.shape
    assert maxShiftDist <= k

    simMat = scipy.sparse.csr_matrix((n, n), dtype=bool)
    for shift in range(1, maxShiftDist + 1):
        simMatUpdate = compare_rows_hash(bpArr[:, shift:], bpArr[:, :-shift])
        simMat = simMat + simMatUpdate
    print("Time until adjacency mat constructed", time.time() - start)

    # from this similarity matrix, generate clusters
    G = nx.from_numpy_array(simMat)
    assemblies = list(nx.connected_components(G))
    print("Time until networkx done", time.time() - start)

    # order clusters by size
    assemblies.sort(key=len, reverse=True)

    return [[anchLst[i] for i in list(cc)] for cc in assemblies]  # output sequences


def main():
    args = parse_args()
    with open(args.input, "r") as f:
        anchLst = [line.strip() for line in f]
    clusters = clusterAnchors(anchLst, args.maxShiftDist)
    # return the first anchor in the cluster as the representative anchor
    with open(args.output, "w") as f:
        for i in range(len(clusters)):
            for j in range(len(clusters[i])):
                f.write(str(i) + "\t" + clusters[i][j] + "\n")
    return None


if __name__ == "__main__":
    main()
