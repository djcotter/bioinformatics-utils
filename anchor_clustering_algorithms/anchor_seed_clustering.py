"""
anchor_seed_clustering.py

Take in a list of anchor sequences, identify the most unique N anchors and create a seed dictionary,
then assign each anchor to the seed that it is most similar to, providing a ranking of the most
similar anchors to each seed. Output a two column file with the first column being the seed id and the
second column being the anchor sequence.

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
        description="Cluster anchors based off of similarity metric."
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
        "-k",
        "--num_clusters",
        type=int,
        default=5000,
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
    return anchors


# create a seed dictionary from the most unique N anchors
def create_seed_dict(anchors, N, distance_threshold=10):
    # create a list to store the seeds
    unique_seeds = []
    # place the first anchor in the seed dictionary
    unique_seeds.append(anchors[0])
    # proceed to the next anchor and compare it to the seeds
    # if it is unique enough, add it to the seed dictionary
    # first shuffle the anchors to avoid bias
    np.random.shuffle(anchors)
    for i in range(1, len(anchors)):
        unique = True
        for seed in unique_seeds:
            if Levenshtein.distance(anchors[i], seed) < distance_threshold:
                unique = False
                break
        if unique:
            unique_seeds.append(anchors[i])
        if len(unique_seeds) >= N:
            break
    # create a seed dictionary
    seed_dict = {i: seed for i, seed in enumerate(unique_seeds)}
    return seed_dict


# for each seed, calculate each anchors similarity and return a ranked
# list of the most similar anchors
def assign_anchors_to_seeds(anchors, seed_dict, distance_threshold=9):
    # create a dictionary to store the cluster assignments
    rankings = {i: [] for i in seed_dict.keys()}
    # loop through the anchors and calculate the similarity to each seed
    for i, anchor in enumerate(anchors):
        best_seed = None
        best_distance = distance_threshold
        for seed_id, seed in seed_dict.items():
            distance = Levenshtein.distance(anchor, seed)
            if distance < best_distance:
                best_seed = seed_id
                best_distance = distance
        if best_seed is not None:
            rankings[best_seed].append(anchor)
        if i % 1000 == 0:
            print(f"Processed {i} anchors")
    return rankings


# main function
def main():
    # parse arguments
    args = parse_args()
    # read in the anchor sequences from the input file
    anchors = read_anchors(args.input_file)
    # create a seed dictionary from the most unique N anchors
    seed_dict = create_seed_dict(anchors, args.num_clusters, args.distance_threshold)
    # assign each anchor to the seed that it is most similar to
    rankings = assign_anchors_to_seeds(anchors, seed_dict, args.distance_threshold - 1)
    # write the cluster assignments to the output file
    with open(args.output_file, "w") as f:
        for seed_id, top_anchors in rankings.items():
            for anchor in top_anchors:
                # write the seed id and the anchor sequence
                # print(f"{seed_id}\t{anchor}")
                f.write(f"{seed_id}\t{anchor}")


if __name__ == "__main__":
    main()
