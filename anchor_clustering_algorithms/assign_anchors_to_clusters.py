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
from functools import partial
import multiprocessing
from multiprocessing import Pool
import time


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
    parser.add_argument(
        "--num_cpus",
        type=int,
        default=None,
        help="Number of CPUs to use for parallel processing.",
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
    anchors = [anchor.strip() for anchor in anchors]
    return anchors


# define a function to process each anchor and assign it to a cluster
def process_anchor(
    anchor,
    cluster_assignments,
    metric,
    distance_threshold,
):
    # initialize dictionary to store the similarity scores for each cluster
    similarity_scores = {}
    # loop through each cluster and calculate the similarity score
    for cluster_id, assigned_anchors in cluster_assignments.items():
        # randomly select up to 5 anchors from the cluster to compare to
        assigned_anchors = np.random.choice(
            assigned_anchors, min(5, len(assigned_anchors))
        )
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
        return None
    assigned_cluster = min(similarity_scores, key=lambda x: similarity_scores[x])
    return assigned_cluster


# define a helper function to print the progress
def print_progress(results, total_anchors, start_time):
    num_anchors_finished = len(results)
    time_passed = time.time() - start_time
    progress = f"Anchors finished: {num_anchors_finished}/{total_anchors}. Time passed: {time_passed:.2f} seconds."
    print(progress, end="\r", flush=True)


# assign new anchor sequences to clusters based on similarity metric
def assign_anchors_to_clusters(
    cluster_assignments,
    anchors,
    metric="lev",
    distance_threshold=9,
    num_cpus=None,
):  # loop through each anchor and calculate the similarity to each cluster

    # parallelize the loop using multiprocessing.Pool
    process_anchor_partial = partial(
        process_anchor,
        cluster_assignments=cluster_assignments,
        metric=metric,
        distance_threshold=distance_threshold,
    )

    # detect the number of available CPUs
    if num_cpus is not None:
        print(f"Using {num_cpus} CPUs.")
    else:
        try:
            num_cpus = multiprocessing.cpu_count()
            print(f"Detected {num_cpus} CPUs.")
        except Exception as e:
            print(e)
            print("Could not detect number of CPUs. Defaulting to 1 CPU.")
            num_cpus = 1

    # create a pool of workers to process the anchors
    with Pool(num_cpus) as pool:
        start_time = time.time()
        results = []
        for i, result in enumerate(
            pool.imap_unordered(process_anchor_partial, anchors, chunksize=5000)
        ):
            results.append(result)
            print_progress(results, len(anchors), start_time)
        print()  # print a new line after the progress bar

    # create a dictionary to store the cluster assignments
    rankings = {cluster_id: [] for cluster_id in cluster_assignments.keys()}
    # loop through the results and assign the anchors to the clusters
    for i, assigned_cluster in enumerate(results):
        if assigned_cluster is not None:
            rankings[assigned_cluster].append(anchors[i])
    return rankings


# write the cluster assignments to the output file
def write_cluster_assignments(output_file, rankings):
    with open(output_file, "w") as f:
        for cluster_id, assigned_anchors in rankings.items():
            for anchor in assigned_anchors:
                f.write(f"{cluster_id}\t{anchor}\n")


# main function
def main():
    args = parse_args()
    print("Reading anchor sequences and cluster assignments...")
    cluster_assignments = read_anchors(args.input_file)
    new_anchors = read_new_anchors(args.new_anchors)
    print("Assigning new anchor sequences to clusters...")
    rankings = assign_anchors_to_clusters(
        cluster_assignments,
        new_anchors,
        metric=args.metric,
        distance_threshold=args.distance_threshold,
        num_cpus=args.num_cpus,
    )
    print("Writing cluster assignments to output file...")
    write_cluster_assignments(args.output_file, rankings)


if __name__ == "__main__":
    main()
