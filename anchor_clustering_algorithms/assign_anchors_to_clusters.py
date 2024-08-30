"""
assign_anchors_to_clusters.py

Given a two column file with the first column being the cluster id and the second column being the anchor sequence,
and a new list of anchor sequences, assign the new anchor sequences to the clusters based off of a chosen similarity metric,
calculating the similarity between the new anchor sequences and a representative from each cluster.

Daniel Cotter
08/29/2024
"""

# Import necessary libraries
import numpy as np
import Levenshtein
import argparse


# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign anchors to clusters based off of similarity metric."
    )
    parser.add_argument(
        "-m",
        "--metric",
        type=str,
        default="lev",
        help="Similarity metric to use for clustering.",
    )
    parser.add_argument(
        "-t",
        "--distance_threshold",
        type=int,
        default=10,
        help="Threshold for similarity metric.",
    )
    # add positional argument for input file
    parser.add_argument(
        "input_file",
        type=str,
        help="Input file containing the list of anchor sequences and their cluster assignments.",
    )
    # add positional argument for new anchor sequences
    parser.add_argument(
        "new_anchors",
        type=str,
        help="File containing the list of new anchor sequences to assign to clusters.",
    )
    # add positional argument for output file
    parser.add_argument(
        "output_file",
        type=str,
        help="Output file containing the cluster assignments for the new anchor sequences.",
    )
    return parser.parse_args()


# read in the anchor sequences and their cluster assignments from the input file
def read_anchors(input_file):
    with open(input_file, "r") as f:
        lines = f.readlines()
    cluster_assignments = {}
    for line in lines:
        cluster_id, anchor = line.strip().split("\t")
        if cluster_id in cluster_assignments:
            cluster_assignments[cluster_id].append(anchor)
        else:
            cluster_assignments[cluster_id] = [anchor]
    return cluster_assignments


# read in the new anchor sequences from the input file
def read_new_anchors(new_anchors):
    with open(new_anchors, "r") as f:
        anchors = f.readlines()
    return anchors


# assign new anchor sequences to clusters based on similarity metric
def assign_anchors_to_clusters(
    cluster_assignments,
    anchors,
    metric="lev",
    distance_threshold=9,
    anchors_per_cluster=10,
):
    rankings = {}
    # loop through the old clusters and calculate the similarity between the new anchors and the old cluster representatives
    for cluster_id, old_anchors in cluster_assignments.items():
        cluster_representative = old_anchors[0]
        similarities = []
        for anchor in anchors:
            if metric == "lev":
                dist = Levenshtein.distance(
                    cluster_representative, anchor, score_cutoff=distance_threshold + 1
                )
            else:
                raise ValueError("Invalid metric specified.")
            if dist < distance_threshold:
                similarities.append(dist)
        # sort the similarities and store the top N
        top_anchors = np.argsort(similarities)[:anchors_per_cluster]
        rankings[cluster_id] = [anchors[i] for i in top_anchors]
    return rankings


# write the cluster assignments to the output file
def write_cluster_assignments(output_file, rankings):
    with open(output_file, "w") as f:
        for cluster_id, assigned_anchors in rankings.items():
            for anchor in assigned_anchors:
                f.write(f"{cluster_id}\t{anchor}")


# main function
def main():
    args = parse_args()
    cluster_assignments = read_anchors(args.input_file)
    new_anchors = read_new_anchors(args.new_anchors)
    rankings = assign_anchors_to_clusters(
        cluster_assignments,
        new_anchors,
        metric=args.metric,
        distance_threshold=args.distance_threshold,
        anchors_per_cluster=10,
    )
    write_cluster_assignments(args.output_file, rankings)


if __name__ == "__main__":
    main()
