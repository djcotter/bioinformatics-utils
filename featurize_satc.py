"""
featurize_satc.py
Daniel Cotter
04/17/2024

This script takes in a file of sample, anchor, target, count
as well as a list of anchors to care about and outputs a fasta
file where the sample is the header and the sequence is composed
of the target for each anchor in the anchor file in order. If a
target is not present, the sequence is composed of Ns.
"""

# import statements
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# function definitions
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script takes in a file of sample, anchor, target, count as well as a list of anchors to care about"
        "and outputs a fasta file where the sample is the header and the sequence is composed of the target for each anchor in the anchor file in order."
        "If a target is not present, the sequence is composed of Ns."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="The input file containing the sample, anchor, target, count",
    )
    parser.add_argument(
        "anchor_file", type=str, help="The file containing the anchors in order"
    )
    parser.add_argument("output_file", type=str, help="The output fasta file")
    parser.add_argument(
        "-a",
        "--include_anchor",
        action="store_true",
        help="Include the anchor in the output fasta file",
    )
    return parser.parse_args()


def read_anchors(anchor_file):
    """
    Reads in the anchor file and returns a list of the anchors
    """
    anchors = []
    with open(anchor_file, "r") as f:
        for line in f:
            anchors.append(line.strip())
    return anchors


def prefilter_input(input_file, anchor_file, tmp_file):
    """
    Uses grep to filter the input file to only include the
    anchors. Then sorts the input file by sample.
    Saves the output to a temporary file.
    """
    # filter the input file to only include the anchors
    os.system(f"grep -Ff {anchor_file} {input_file} > {tmp_file}")
    # sort the input file by sample
    os.system(f"sort -k1,1 -o {tmp_file} {tmp_file}")
    return None


def write_fasta(
    input_file, anchors, output_file, missing_token="NNNNN", include_anchor=False
):
    """
    Reads in the input file and writes the fasta file in the anchor order
    of the anchor file. If a target is not present, the sequence is composed
    of Ns.
    """
    with open(output_file, "w") as f:
        with open(input_file, "r") as g:
            sample_count = 0
            current_sample = ""
            current_anchors = {}
            for line in g:
                sample, anchor, target, count = line.strip().split("\t")[:4]
                if sample[1] == "anchor":
                    continue
                if sample != current_sample:
                    if current_sample != "":
                        sample_count += 1
                        if sample_count % 1000 == 0:
                            print(f"Processed {sample_count} samples")
                        # write the current sample to the fasta file
                        if include_anchor:
                            # make the fasta format anchor followed by target
                            seq = "".join(
                                [
                                    anchor
                                    + current_anchors.get(anchor, missing_token)[0]
                                    for anchor in anchors
                                ]
                            )
                        else:
                            seq = "".join(
                                [
                                    current_anchors.get(anchor, missing_token)[0]
                                    for anchor in anchors
                                ]
                            )
                        record = SeqRecord(Seq(seq), id=current_sample, description="")
                        SeqIO.write(record, f, "fasta")
                    # update the current sample
                    current_sample = sample
                    current_anchors = {}
                # update the current anchors by comparing the counts
                if anchor not in current_anchors:
                    current_anchors[anchor] = [target, count]
                elif count > current_anchors[anchor][1]:
                    current_anchors[anchor] = [target, count]
            # write the last sample to the fasta file
            if include_anchor:
                # make the fasta format anchor followed by target
                seq = "".join(
                    [
                        anchor + current_anchors.get(anchor, missing_token)[0]
                        for anchor in anchors
                    ]
                )
            else:
                seq = "".join(
                    [
                        current_anchors.get(anchor, missing_token)[0]
                        for anchor in anchors
                    ]
                )
            record = SeqRecord(Seq(seq), id=current_sample, description="")
            SeqIO.write(record, f, "fasta")
    return None


def main():
    # parse arguments
    args = parse_args()
    input_file = args.input_file
    anchor_file = args.anchor_file
    output_file = args.output_file
    include_anchor = args.include_anchor
    # read in the anchors
    anchors = read_anchors(anchor_file)
    # filter the input file
    # temp file should be basename input file with .tmp extension
    tmp_file = f"{input_file}.tmp"
    prefilter_input(input_file, anchor_file, tmp_file)
    # write the fasta file
    write_fasta(tmp_file, anchors, output_file, include_anchor=include_anchor)
    # remove the temporary file
    os.system(f"rm {tmp_file}")
    return None


# execute main
if __name__ == "__main__":
    main()
