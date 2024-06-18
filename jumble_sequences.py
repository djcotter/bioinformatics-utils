"""
jumble_sequences.py

This script takes in a fasta file composed of sequences to be embedded in a DNA LLM
it then performs an operation (from a set of operations) to jumble/mutate the sequences 
and writes them out to a new fasta file. The operations are as follows:
    - Remove any Ns from the sequence (--remove_n)
    - reverse the entire sequence (--reverse)
    - truncate the sequence to a specified length (--truncate <int>)
    - mutate the sequence by changing all A's to T's (--mutate_a_to_t)
    - mutate the sequence by changing all C's to G's (--mutate_c_to_g)
    - mutate the sequence by changing a specified fraction of the bases (--mutate <float>)

Also takes an argument to output as a tsv file instead of a fasta file with
the sequence name in the first column and the sequence in the second column.

Usage:
python jumble_sequences.py input.fasta output.fasta --remove_n --reverse --truncate 100
--mutate_a_to_t --mutate_c_to_g --mutate 0.1 --output_tsv INPUT_FILE OUTPUT_FILE


Daniel Cotter - 06/18/2024
"""

# Import Statements
from Bio import SeqIO
import argparse
import random
import gzip
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Function Definitions
def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(description="Jumble sequences in a fasta file")
    parser.add_argument("input", help="Input fasta file")
    parser.add_argument("output", help="Output file")
    parser.add_argument(
        "--remove_n", help="Remove Ns from the sequence", action="store_true"
    )
    parser.add_argument("--reverse", help="Reverse the sequence", action="store_true")
    parser.add_argument(
        "--truncate",
        help="Truncate the sequence to a specified length",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--mutate_a_to_t", help="Mutate all A's to T's", action="store_true"
    )
    parser.add_argument(
        "--mutate_c_to_g", help="Mutate all C's to G's", action="store_true"
    )
    parser.add_argument(
        "--mutate_fraction",
        help="Mutate a specified fraction of the bases",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--output_tsv",
        help="Output as a tsv file instead of a fasta file",
        action="store_true",
    )
    return parser.parse_args()


def import_fasta_records(input_file: str) -> list:
    """
    Import all fasta records as a list from a fasta file
    """
    if "fastq" in input_file or "fq" in input_file:
        filetype = "fastq"
    elif "fasta" in input_file or "fa" in input_file:
        filetype = "fasta"
    try:
        with gzip.open(input_file, "rt") as f:
            records = list(SeqIO.parse(f, format=filetype))
    except gzip.BadGzipFile:
        with open(input_file, "rt") as f:
            records = list(SeqIO.parse(f, format=filetype))
    return records


def remove_n(record: str) -> str:
    """
    Remove Ns from the sequence
    """
    return str(record).replace("N", "").replace("n", "")


def reverse(record: str) -> str:
    """
    Reverse the sequence
    """
    return str(record)[::-1]


def truncate(record: str, length: int) -> str:
    """
    Truncate the sequence to a specified length
    """
    return str(record)[:length]


def mutate_a_to_t(record: str) -> str:
    """
    Mutate all A's to T's
    """
    return str(record).replace("A", "T").replace("a", "t")


def mutate_c_to_g(record: str) -> str:
    """
    Mutate all C's to G's
    """
    return str(record).replace("C", "G").replace("c", "g")


def mutate_fraction(record: str, fraction: float) -> str:
    """
    Mutate a specified fraction of the bases
    Mutations are currently made A->T, T->A, C->G, G->C
    """
    record = list(record)
    for i in range(len(record)):
        if random.random() < fraction:
            if record[i] == "A":
                record[i] = "T"
            elif record[i] == "T":
                record[i] = "A"
            elif record[i] == "C":
                record[i] = "G"
            elif record[i] == "G":
                record[i] = "C"
    return "".join(record)


def main() -> None:
    """
    Main function
    """
    # parse arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output
    records = import_fasta_records(input_file)
    print(records)
    jumbled_records = []
    for record in records:
        name = record.id
        record = record.seq
        if args.remove_n:
            record = remove_n(record)
        if args.reverse:
            record = reverse(record)
        if args.truncate:
            record = truncate(record, args.truncate)
        if args.mutate_a_to_t:
            record = mutate_a_to_t(record)
        if args.mutate_c_to_g:
            record = mutate_c_to_g(record)
        if args.mutate_fraction:
            record = mutate_fraction(record, args.mutate_fraction)
        jumbled_records.append(SeqRecord(Seq(record), id=name, description=""))
    if args.output_tsv:
        with open(output_file, "w") as f:
            csv_writer = csv.writer(f, delimiter="\t")
            csv_writer.writerows(
                [[record.id, record.seq] for record in jumbled_records]
            )
    else:
        with open(output_file, "w") as f:
            SeqIO.write(jumbled_records, f, "fasta")
    return None


# Main
if __name__ == "__main__":
    main()
