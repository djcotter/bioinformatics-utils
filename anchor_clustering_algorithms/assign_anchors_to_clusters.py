"""
assign_anchors_to_clusters.py

Given a two column file with the first column being the cluster id and the second column being the anchor sequence,
and a new list of anchor sequences, assign the new anchor sequences to the clusters based off of a chosen similarity metric,
calculating the similarity between the new anchor sequences and a representative from each cluster.

Daniel Cotter
08/29/2024
"""

# Import necessary libraries
# import numpy as np
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
):
    # initialize dictionary to store the rankings of the new anchors for each cluster
    rankings = {}
    count = 0
    # loop through each anchor and calculate the similarity to each cluster
    for anchor in anchors:
        # initialize dictionary to store the similarity scores for each cluster
        similarity_scores = {}
        # loop through each cluster and calculate the similarity score
        for cluster_id, assigned_anchors in cluster_assignments.items():
            # randomly select up to 5 anchors from the cluster to compare to
            # assigned_anchors = np.random.choice(
            #     assigned_anchors, min(5, len(assigned_anchors))
            # )
            assigned_anchors = assigned_anchors[0 : min(5, len(assigned_anchors))]
            # calculate the similarity score based on the chosen metric
            if metric == "lev":
                similarity = min(
                    [
                        Levenshtein.distance(
                            anchor, assigned_anchor, score_cutoff=distance_threshold + 1
                        )
                        for assigned_anchor in assigned_anchors
                    ]
                )
            else:
                raise ValueError("Invalid metric specified.")
            # store the similarity score in the dictionary
            similarity_scores[cluster_id] = similarity
        # remove any clusters that are not similar enough
        similarity_scores = {
            cluster_id: similarity
            for cluster_id, similarity in similarity_scores.items()
            if similarity <= distance_threshold
        }
        # assign the anchor to the cluster with the highest similarity score if there is one
        # otherwise don't include the anchor in any cluster
        if similarity_scores == {}:
            continue
        assigned_cluster = min(similarity_scores, key=lambda x: similarity_scores[x])
        # add the anchor to the rankings dictionary
        if assigned_cluster in rankings:
            rankings[assigned_cluster].append(anchor)
        else:
            rankings[assigned_cluster] = [anchor]
        # print progress
        count = count + 1
        if count % 1000 == 0:
            print(f"Processed {count} anchors")
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
    )
    write_cluster_assignments(args.output_file, rankings)


if __name__ == "__main__":
    main()
