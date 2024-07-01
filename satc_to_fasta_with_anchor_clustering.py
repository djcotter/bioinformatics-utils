"""
satc_to_fasta_with_anchor_clustering.py

Convert a satc file to a fasta file with anchor clustering. SATC files are 
in the format of: Sample, Anchor, Target, Count.

This script will take in a SATC file and a two column anchor file where column
1 is the cluster id and column 2 is the anchor. It will then convert the SATC file
into a fasta file picking one representative anchor from each cluster. For example, if sample 2
does not have the first anchor in cluster 1 but does have the second anchor in cluster 1, the second
anchor will be chosen as the representative anchor for that cluster for sample 2.

Usage:
python satc_to_fasta_with_anchor_clustering.py [-n INT] [-d anchor_file] [-i satc_dump_file] [-o output_file]

Daniel Cotter - 07/01/2024
"""

# Import Statements
from Bio import SeqIO
import argparse
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Function Definitions
def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(
        description="Convert a satc file to a fasta file with anchor clustering"
    )
    parser.add_argument(
        "-n",
        "--num_targets",
        help="Number of targets to ouput per anchor",
        default=1,
        type=int,
    )
    parser.add_argument("-d", "--anchor_file", help="Anchor file", required=True)
    parser.add_argument("-i", "--satc_dump_file", help="SATC dump file", required=True)
    parser.add_argument("-o", "--output_file", help="Output fasta file", required=True)
    return parser.parse_args()


def read_anchors(anchor_file):
    """
    Reads in the anchor file and returns a list of the anchors
    """
    anchors_dict = {}
    with open(anchor_file, "r") as f:
        for line in f:
            cluster, anchor = line.strip().split()
            anchors_dict[anchor] = cluster
    return anchors_dict


def prefilter_input(input_file, anchor_file, tmp_file):
    """
    Uses grep to filter the input file to only include the
    anchors. Then sorts the input file by sample.
    Saves the output to a temporary file.
    """
    # filter the input file to only include the anchors
    os.system(f"grep -Ff <(cut -f2 {anchor_file}) {input_file} > {tmp_file}")
    # sort the input file by sample
    os.system(f"sort -k1,1 -o {tmp_file} {tmp_file}")
    return None


def write_fasta(input_file, anchor_dict, output_file, num_targets=1, anchor_length=27):
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
            anchor_clusters = []
            for line in g:
                sample, anchor, target, count = line.strip().split("\t")[:4]
                count = int(count)
                if anchor == "anchor":
                    continue
                if sample != current_sample:
                    if current_sample != "":
                        sample_count += 1
                        if sample_count % 1000 == 0:
                            print(f"Processed {sample_count} samples")
                        # write the current sample to the fasta file
                        seq = ""
                        for anchor in current_anchors:
                            new_seq = "".join(current_anchors[anchor][:num_targets])
                            if len(current_anchors[anchor]) < num_targets:
                                new_seq += (
                                    "N"
                                    * anchor_length
                                    * (num_targets - len(current_anchors[anchor]))
                                )
                            seq += new_seq
                        record = SeqRecord(Seq(seq), id=current_sample, description="")
                        SeqIO.write(record, f, "fasta")
                    current_sample = sample
                    current_anchors = []
                    anchor_clusters = []
                if anchor_dict[anchor] not in anchor_clusters:
                    anchor_clusters.append(anchor_dict[anchor])
                    current_anchors[anchor] = [target]
                elif anchor in current_anchors:
                    if len(current_anchors[anchor]) < num_targets:
                        current_anchors[anchor].append(target)
            # write the last sample to the fasta file
            seq = ""
            for anchor in current_anchors:
                new_seq = "".join(current_anchors[anchor][:num_targets])
                if len(current_anchors[anchor]) < num_targets:
                    new_seq += (
                        "N"
                        * anchor_length
                        * (num_targets - len(current_anchors[anchor]))
                    )
                seq += new_seq
            record = SeqRecord(Seq(seq), id=current_sample, description="")
            SeqIO.write(record, f, "fasta")
    return None


# Main
def main():
    """
    Main function
    """
    # parse arguments
    args = parse_args()
    num_targets = args.num_targets
    anchor_file = args.anchor_file
    satc_dump_file = args.satc_dump_file
    output_file = args.output_file
    # read in the anchor file
    anchor_dict = read_anchors(anchor_file)
    # prefilter the input file
    tmp_file = "anchors_tmp.txt"
    prefilter_input(satc_dump_file, anchor_file, tmp_file)
    # write the fasta file
    write_fasta(tmp_file, anchor_dict, output_file, num_targets)
    # remove the temporary file
    os.system(f"rm {tmp_file}")
    return None


# execute main
if __name__ == "__main__":
    main()
