"""
cluster_anchors.py

Take in a list of anchor sequences and cluster them based off of a provided similarity metric.
Build an adjecency matrix based off of a chosen similarity metric and then use spectral 
clustering to cluster the anchors. Output the a two column file with the first column being the
cluster id and the second column being the anchor sequence.

Daniel Cotter
08/27/2024
"""

# Import necessary libraries
import time
import numpy as np
import scipy
import nltk
import argparse
from sklearn.cluster import SpectralClustering


# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster anchors based off of similarity metric."
    )
    parser.add_argument(
        "-m",
        "--metric",
        type=str,
        required=True,
        default="lev",
        help="Similarity metric to use for clustering.",
    )
    parser.add_argument(
        "-t",
        "--distance_threshold",
        type=int,
        default=0,
        help="Threshold for similarity metric.",
    )
    parser.add_argument(
        "-k",
        "--num_clusters",
        type=int,
        required=True,
        help="Number of clusters to partition the anchors into.",
    )
    # add positional argument for input file
    parser.add_argument(
        "input_file",
        type=str,
        help="Input file containing the list of anchor sequences.",
    )
    # add positional argument for output file
    parser.add_argument(
        "output_file",
        type=str,
        help="Output file containing the cluster assignments.",
    )
    return parser.parse_args()


# read in the anchor sequences from the input file
def read_anchors(input_file):
    with open(input_file, "r") as f:
        anchors = f.readlines()
    return [anchor.strip() for anchor in anchors]


# calculate an adjacency matrix based off of the chosen similarity metric
def adjacency_matrix(anchors, metric, threshold=5):
    start_time = time.time()
    if metric == "lev":
        simMat = scipy.sparse.csr_matrix((len(anchors), len(anchors)), dtype=np.uint8)
        for i in range(len(anchors)):
            for j in range(i, len(anchors)):
                # levenshtien distance
                dist = nltk.edit_distance(anchors[i], anchors[j])
                if dist <= threshold:
                    simMat[i, j] = 1
                    simMat[j, i] = 1
        print(
            "Time until adjacency mat constructed: {}".format(time.time() - start_time)
        )
        return simMat
    else:
        raise ValueError("Invalid similarity metric. Use 'lev'.")


# cluster the anchors using spectral clustering
def cluster_anchors(simMat, num_clusters):
    start = time.time()
    sc = SpectralClustering(num_clusters, affinity="precomputed")
    cluster_assignments = sc.fit_predict(simMat)
    print("Time until clustering done: {}".format(time.time() - start))
    return cluster_assignments


# write the cluster assignments to the output file
def write_clusters(output_file, cluster_assignments, anchors):
    with open(output_file, "w") as f:
        for i, cluster_id in enumerate(cluster_assignments):
            f.write(f"{cluster_id}\t{anchors[i]}\n")


# main function
def main():
    args = parse_args()
    anchors = read_anchors(args.input_file)
    if args.distance_threshold > 0:
        simMat = adjacency_matrix(anchors, args.metric, args.distance_threshold)
    else:
        # default to 5 (assuming lev here)
        simMat = adjacency_matrix(anchors, args.metric, 5)
    cluster_assignments = cluster_anchors(simMat, args.num_clusters)
    print(cluster_assignments)
    # write_clusters(args.output_file, cluster_assignments, anchors)


if __name__ == "__main__":
    main()
