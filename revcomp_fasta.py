"""
revcomp_fasta.py
Daniel Cotter -- 04/05/24

Get the reverse complement of all records in a fasta file and write them to a new file
"""

# Import Statements
from Bio import SeqIO
import argparse


def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(description="Get reverse complement of fasta file")
    parser.add_argument("input", help="Input fasta file")
    parser.add_argument("output", help="Output fasta file")
    return parser.parse_args()


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    with open(output_file, "w") as f:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.reverse_complement()
            SeqIO.write(record, f, "fasta")


# run main function
if __name__ == "__main__":
    main()
