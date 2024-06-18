"""
kmer_decomposition.py

Get all kmers of a defined length (-k) by sliding in windows (with offset -s; default is no overlap)
along all records in a fasta file. Record only unique kmers and keep track of where they occur in the 
original sequences. Write a fasta composed of only the unique kmers and a tsv file composed of kmer ids
that can be used to map the kmers back to the original sequences.

Daniel Cotter -- 06/18/24
"""

# Import Statements
from Bio import SeqIO
import argparse
from more_itertools import windowed
import gzip
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Function Definitions
def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(description="Get kmers from a fasta file")
    parser.add_argument("input", help="Input fasta file")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("-k", "--kmer", help="Kmer length", default=27, type=int)
    parser.add_argument("-s", "--step", help="Step size", default=None, type=int)
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


def get_kmers(record: str, k: int, s: int) -> list:
    """
    Return list of all kmers of length k with step size s from current record
    Will not report the last kmer if the step size causes the window to advance
    past the end of the record
    """
    # list of kmers with None values removed
    kmers = windowed(str(record), k, step=s)
    kmers = ["".join(kmer) for kmer in kmers if None not in kmer]
    # list of indices in the sequence where each kmer starts
    [start, end] = zip(*[(i, i + k) for i in range(0, len(record) - k + 1, s)])
    return kmers, start, end


def get_unique_kmers(current_kmer_dict: dict, new_kmers: list) -> dict:
    """
    Take in an existing dictionary of kmers and their unique indices and add new kmers
    with a new index if thet are not already in the dictionary
    """
    new_i = len(current_kmer_dict)
    for kmer in new_kmers:
        if kmer not in current_kmer_dict:
            current_kmer_dict[kmer] = new_i
            new_i += 1
    return current_kmer_dict


def generate_kmer_order(kmers: list, start: list, end: list, kmer_dict: dict) -> list:
    """
    Generate a list of kmer ids that correspond to the kmers in the original sequence
    """
    kmer_order = [kmer_dict[kmer] for kmer in kmers]
    return zip(kmers, kmer_order, start, end)


def write_out_fasta(sequences: dict, output_file: str) -> None:
    """
    Write out fasta file with kmer sequences
    """
    with open(output_file + "_unique_kmers.fasta", "w") as f:
        for kmer in sequences:
            record = SeqRecord(Seq(kmer), id=str(sequences[kmer]), description="")
            SeqIO.write(record, f, "fasta")


def write_out_tsv(sequences: list, output_file: str) -> None:
    """
    Write out tsv file with 4 columns: record id, kmer id, start, end
    """
    with open(output_file + "_kmer_ordering.tsv", "w") as f:
        csv_writer = csv.writer(f, delimiter="\t")
        csv_writer.writerows(sequences)


def main() -> None:
    """
    Main function
    """
    # parse arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output_prefix
    kmer_length = args.kmer
    if args.step:
        step_size = args.step
    else:
        step_size = kmer_length
    # import fasta records
    records = import_fasta_records(input_file)
    # get kmers and indices
    kmer_dict = {}
    kmer_order = []
    for record in records:
        kmers, start, end = get_kmers(record.seq, kmer_length, step_size)
        kmer_dict = get_unique_kmers(kmer_dict, kmers)
        kmer_ordering = generate_kmer_order(kmers, start, end, kmer_dict)
        kmer_ordering = [list(km) for km in kmer_ordering]
        kmer_ordering = [[record.id] + km for km in kmer_ordering]
        # add the record id to each entry in the kmer_order list
        kmer_order.extend(kmer_ordering)
    # write out fasta and tsv
    output_file = output_file + f"_k{str(kmer_length)}_s{str(step_size)}"
    write_out_fasta(kmer_dict, output_file)
    write_out_tsv(kmer_order, output_file)
    return None


if __name__ == "__main__":
    main()
