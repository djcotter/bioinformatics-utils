"""
fasta_to_tsv.py

Convert a fasta file to a tsv file with 2 columns: record id, sequence
or convert a tsv file with 2 columns: record id, sequence to a fasta file

Usage:
python fasta_to_tsv.py input output

Detect input file format by extension: .fasta, .fa, .fastq, .fq, .tsv
Detect output file format by extension: .fasta, .fa, .fastq, .fq, .tsv

Daniel Cotter 06/18/2024
"""

# Import Statements
from Bio import SeqIO
import argparse
import gzip
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Function Definitions
def parse_args():
    """
    Parse commandline arguments
    """
    parser = argparse.ArgumentParser(description="Convert fasta to tsv or tsv to fasta")
    parser.add_argument("input", help="Input file")
    parser.add_argument("output", help="Output file")
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


def write_fasta(records: list, output_file: str):
    """
    Write all records to a fasta file
    """
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")


def write_tsv(records: list, output_file: str):
    """
    Write all records to a tsv file
    """
    with open(output_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for record in records:
            writer.writerow([record.id, str(record.seq)])


def import_tsv_records(input_file: str) -> list:
    """
    Import all records as a list from a tsv file
    """
    with open(input_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        records = []
        for row in reader:
            record = SeqRecord(Seq(row[1]), id=row[0], description="")
            records.append(record)
    return records


# Main
def main():
    args = parse_args()
    if (
        "fasta" in args.input
        or "fa" in args.input
        or "fastq" in args.input
        or "fq" in args.input
    ):
        records = import_fasta_records(args.input)
        if (
            "fasta" in args.output
            or "fa" in args.output
            or "fastq" in args.output
            or "fq" in args.output
        ):
            write_fasta(records, args.output)
        elif "tsv" in args.output:
            write_tsv(records, args.output)
        else:
            print("Output file must be a fasta, fa, fastq, fq, or tsv file")
    elif "tsv" in args.input:
        records = import_tsv_records(args.input)
        if (
            "fasta" in args.output
            or "fa" in args.output
            or "fastq" in args.output
            or "fq" in args.output
        ):
            write_fasta(records, args.output)
        elif "tsv" in args.output:
            write_tsv(records, args.output)
        else:
            print("Output file must be a fasta, fa, fastq, fq, or tsv file")
    else:
        print("Input file must be a fasta, fa, fastq, fq, or tsv file")


if __name__ == "__main__":
    main()
