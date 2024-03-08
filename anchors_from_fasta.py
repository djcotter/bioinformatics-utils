"""
anchors_from_fasta.py

Get anchors of a defined length (-k) by sliding in windows (with offset -s)
along all records (or a single record; --record) in a fasta file and write them to a file
Can also get the reverse complement of the anchors (--reverse)

Daniel Cotter -- 01/16/24
"""

# Import Statements
from Bio import SeqIO
import argparse
from more_itertools import windowed
import gzip
import csv


def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(description="Get anchors from a fasta file")
    parser.add_argument("input", help="Input fasta file")
    parser.add_argument("output", help="Output file")
    parser.add_argument("-k", "--kmer", help="Kmer length", default=27, type=int)
    parser.add_argument("-s", "--step", help="Step size", default=1, type=int)
    parser.add_argument("--record", help="Record to get anchors from", required=False)
    parser.add_argument(
        "-r", "--reverse", help="Reverse complement", action="store_true"
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


def get_kmers(record: str, k: int, s: int) -> list:
    """
    Return list of all kmers of length k with step size s from current record
    Will not report the last kmer if the step size causes the window to advance
    past the end of the record
    """
    kmers = windowed(str(record), k, step=s)
    return ["".join(kmer) for kmer in kmers if None not in kmer]


def write_out(anchors: list, output_file: str) -> None:
    """
    Write out anchors to file
    """
    with open(output_file, "w") as f:
        csv_writer = csv.writer(f, delimiter="\t")
        csv_writer.writerows([[anchor.upper()] for anchor in anchors])


def main() -> None:
    """
    Main function
    """
    # parse arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output
    kmer = args.kmer
    step = args.step
    fasta_record = args.record
    # import fasta records
    records = import_fasta_records(input_file)
    # check if fasta record id was provided and filter records
    if fasta_record:
        records = [record.seq for record in records if record.id == fasta_record]
    else:
        records = [record.seq for record in records]
    if args.reverse:
        reverse_records = [record.reverse_complement() for record in records]
        records.extend(reverse_records)
    # get anchors
    anchors = [get_kmers(record, kmer, step) for record in records]
    # collapse anchors into a single list
    anchors = [anchor for sublist in anchors for anchor in sublist]
    write_out(anchors, output_file)


if __name__ == "__main__":
    main()
