"""
parse_anchor_target_counts.py
Daniel Cotter -- 04/04/24

Parse a fasta file and a list of anchors. Create a dictionary of the anchors and their 
reverse complements. For each fasta record, check for all kmers in the anchor dictionary, and
if an anchor exists, print out the fasta record id, the anchor, and the sequence of length
s following that anchor (the target). Also print a column indicating if it's a reverse complement. 
"""

# import packages
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq


# Parse commandline arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Parse fasta file and list of anchors")
    parser.add_argument("fasta", help="Input fasta file", type=str)
    parser.add_argument("anchors", help="Input anchor file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument(
        "-s", "--target_length", help="Length of target sequence", default=5, type=int
    )
    return parser.parse_args()


# create anchor dictionary
def create_anchor_dict(anchors: str) -> dict:
    anchor_dict = {}
    with open(anchors, "r") as f:
        for line in f:
            # first add the anchor as a key with an empty list to the dictionary
            anchor = line.strip()
            anchor_dict[anchor] = []
            # then add the reverse complemnt of the anchor as a key with an empty list
            anchor_dict[str(Seq(anchor).reverse_complement())] = []
    return anchor_dict


# main function
def main():
    # parse arguments
    args = parse_args()
    anchors = args.anchors
    output_file = args.output
    fasta_file = args.fasta
    target_length = args.target_length
    # create anchor dictionary
    anchor_dict = create_anchor_dict(anchors)
    # get k as the length of the first anchor
    k = len(list(anchor_dict.keys())[0])
    # open output file
    with open(output_file, "w") as o:
        writer = csv.writer(o, delimiter="\t")
        # add a header to the output file
        writer.writerow(["fasta_record", "anchor", "target", "reverse_complement"])
        # iterate through fasta records
        for record in SeqIO.parse(fasta_file, "fasta"):
            # get sequence as a string
            sequence = str(record.seq)
            # iterate through all kmers of length 27
            for i in range(len(sequence) - k):
                kmer = sequence[i : i + k]
                # if the kmer is in the anchor dictionary, add the record id, anchor, and target to the output file
                if kmer in anchor_dict:
                    writer.writerow(
                        [
                            record.description,
                            kmer,
                            sequence[i + k : i + k + target_length],
                            False,
                        ]
                    )
            rev_comp = str(Seq(sequence).reverse_complement())
            # iterate through all kmers of length k in the reverse complement of the sequence
            for i in range(len(rev_comp) - k):
                kmer = rev_comp[i : i + k]
                if kmer in anchor_dict:
                    writer.writerow(
                        [
                            record.description,
                            kmer,
                            rev_comp[i + k : i + k + target_length],
                            True,
                        ]
                    )
    return None


# run main function
if __name__ == "__main__":
    main()
