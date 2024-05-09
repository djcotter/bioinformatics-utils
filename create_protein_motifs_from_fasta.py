"""
create_protein_motifs_from_fasta.py
Daniel Cotter - 05/08/2024

Read in a fasta file and create a two column file where the 
first column is the header and the second column is all of the 
N-width substrings of the sequence.
"""

# import statements
import argparse
from Bio import SeqIO


# function definitions
def parse_args():
    parser = argparse.ArgumentParser(
        description="Read in a fasta file and create a two column file where the "
        "first column is the header and the second column is all of the "
        "N-width substrings of the sequence."
    )
    parser.add_argument(
        "fasta_file",
        type=str,
        help="The input fasta file",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="The output file",
    )
    parser.add_argument(
        "-n",
        "--substring_length",
        type=int,
        default=9,
        help="The length of the substrings to extract from the sequences",
    )
    return parser.parse_args()


def create_substrings(fasta_file, output_file, substring_length):
    with open(output_file, "w") as f_out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.id
            sequence = str(record.seq)
            for i in range(len(sequence) - substring_length + 1):
                substring = sequence[i : i + substring_length]
                f_out.write(f"{header}\t{substring}\n")
    return None


def main():
    args = parse_args()
    create_substrings(args.fasta_file, args.output_file, args.substring_length)


if __name__ == "__main__":
    main()
