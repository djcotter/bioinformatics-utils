"""
run_lookup.py
Daniel Cotter - 05/08/2024

This script takes in a file with a list of anchors (any seqeunce)
as well as a path to a prebuilt lookup table and outputs a 3 column
file with the anchor, lookup query, and lookup stats.
"""

# import statements
import os
import argparse
from pathlib import Path


# function definitions
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script takes in a file with a list of anchors (any seqeunce)"
        "as well as a path to a prebuilt lookup table and outputs a 3 column"
        "file with the anchor, lookup query, and lookup stats."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="The input file containing the list of anchors",
    )
    parser.add_argument(
        "lookup_table",
        type=str,
        help="The path to the prebuilt lookup table",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="The output file",
    )
    parser.add_argument(
        "--splash_bin",
        type=Path,
        default=Path("/oak/stanfor/groups/horence/dcotter1/splash-2.6.1"),
    )
    return parser.parse_args()


def read_anchors(anchor_file):
    # read in the anchors from the input file
    with open(anchor_file, "r") as f:
        anchors = [line.strip() for line in f]
    return anchors


def write_temp_fasta(anchors, temp_fasta):
    # use biopython to generate a fasta file where the header is the anchor
    # and the sequence is the anchor
    with open(temp_fasta, "w") as f:
        for anchor in anchors:
            f.write(f">{anchor}\n{anchor}\n")
    return None


def run_lookup(anchor_fasta, lookup_file, output_file, splash_bin):
    # run the lookup command
    lookup_table = Path(splash_bin) / "lookup_table"
    lookup_cmd = (
        f"{lookup_table} query "
        "--kmer_skip 1 --truncate_paths --stats_fmt with_stats "
        "{Path(anchor_fasta).resolve()} {Path(lookup_file).resolve()} {Path(output_file).resolve}"
    )
    print(f"Running command: {lookup_cmd}")
    os.system(lookup_cmd)
    return None


def read_lookup_output(output_file):
    # read in the output file and return the results
    with open(output_file, "r") as f:
        results = [line.strip().split("\t") for line in f]
    return results


def main():
    # parse arguments
    args = parse_args()
    # read in the anchors from the input file
    anchors = read_anchors(args.input_file)
    # write the anchors to a temporary fasta file
    temp_fasta = Path(args.input_file).with_suffix(".fasta")
    write_temp_fasta(anchors, temp_fasta)
    # run the lookup command
    temp_lookup_output = Path(args.input_file).with_suffix(".lookup")
    run_lookup(temp_fasta, args.lookup_table, temp_lookup_output, args.splash_bin)
    # read in the lookup results
    results = read_lookup_output(temp_lookup_output)
    # remove the temporary files
    os.remove(temp_fasta)
    os.remove(temp_lookup_output)
    # write the lookup results with the anchor added on as the first column
    with open(args.output_file, "w") as f:
        for i in range(len(anchors)):
            f.write(f"{anchors[i]}\t{results[i][0]}\t{results[i][1]}\n")
    return None


if __name__ == "__main__":
    main()
