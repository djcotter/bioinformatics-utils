"""
extract_embeddings_from_torch_objects.py

Simple utility to loop through all .pt files in a directory and extract the embeddings from each file. The embeddings are
then all saved to a singly tsv file with the same name as the directory. The embeddings are saved in the order that they
are read from the .pt files.

Daniel Cotter
09/09/2024
"""

# Import necessary libraries
import os
import torch
import argparse


# Parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract embeddings from PyTorch objects and save to tsv file."
    )
    parser.add_argument(
        "input_dir",
        type=str,
        help="Directory containing PyTorch objects (.pt files) with embeddings.",
    )
    return parser.parse_args()


# Extract embeddings from PyTorch object
def extract_embeddings(pt_file):
    embeddings = torch.load(pt_file)
    embeddings = embeddings["mean_embeddings"][33]
    return embeddings


# Main function
def main():
    args = parse_args()
    embeddings = {}
    for file in os.listdir(args.input_dir):
        name = file.split(".")[0]
        if file.endswith(".pt"):
            embeddings[name] = extract_embeddings(os.path.join(args.input_dir, file))

    output_file = args.input_dir + ".tsv"
    output_file = os.path.join(args.input_dir, output_file)
    with open(output_file, "w") as f:
        for key, value in embeddings.items():
            f.write(key + "\t" + "\t".join([str(x) for x in value]) + "\n")


if __name__ == "__main__":
    main()
