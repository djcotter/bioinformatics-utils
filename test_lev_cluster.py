"""
Deduplicate anchors by clustering them based on their similarity.
Tavor's code
"""

import time
import numpy as np
import scipy
import levenshtein
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
    return parser.parse_args()


# main function
# inputs: anchLst, a list of sequences; maxShiftDist, the maximum number of bases that sequences can be shifted and considered in the same cluster
# output: a list of lists, where each list is a cluster of sequences
def clusterAnchors(anchLst):
    start = time.time()

    # construct a sparse scipy matrix using levenshtein distance between anchors
    simMat = scipy.sparse.csr_matrix((len(anchLst), len(anchLst)), dtype=np.uint8)
    for i in range(len(anchLst)):
        for j in range(i, len(anchLst)):
            # levenshtien distance
            dist = levenshtein.distance(anchLst[i], anchLst[j], score_cutoff=7)
            if dist <= 5:
                simMat[i, j] = 1
                simMat[j, i] = 1
    print("Time until adjacency mat constructed", time.time() - start)

    # from this similarity matrix, generate clusters
    start = time.time()
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
    clusters = clusterAnchors(anchLst)
    # return the first anchor in the cluster as the representative anchor
    with open(args.output, "w") as f:
        for i in range(len(clusters)):
            for j in range(len(clusters[i])):
                f.write(str(i) + "\t" + clusters[i][j] + "\n")
    return None


if __name__ == "__main__":
    main()
