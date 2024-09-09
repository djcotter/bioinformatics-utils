"""
translate_fasta.py

Takes in a FASTA file and translates the nucleotide sequences to amino acid sequences. Outputs a new FASTA file with
six reading frames for each nucleotide sequence. The output file will have the same header as the input file with the
reading frame appended to the end of the header. The reading frame will be indicated by a number from 1 to 6.

Daniel Cotter
09/06/2024
"""

# Import necessary libraries
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Parse arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences in a FASTA file to amino acid sequences."
    )
    parser.add_argument(
        "-t",
        "--translation_table",
        type=int,
        default=1,
        help="Translation table to use for translating nucleotide sequences to amino acid sequences.",
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Input FASTA file containing nucleotide sequences to be translated.",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Output FASTA file containing translated amino acid sequences.",
    )
    return parser.parse_args()


# define functions
def translate_sequence(nucleotide_seq, translation_table):
    """
    Translates a nucleotide sequence to an amino acid sequence
    using the specified translation table
    """
    seq = Seq(nucleotide_seq)
    return seq.translate(table=translation_table)


# def translate_fasta(input_file, output_file, translation_table):
#     """
#     Translates nucleotide sequences in a FASTA file to amino acid sequences
#     and writes all 6 reading frames to the output file.
#     """
#     with open(output_file, "w") as f_out:
#         for record in SeqIO.parse(input_file, "fasta"):
#             nucleotide_seq = Seq(str(record.seq))
#             for i in range(3):
#                 frame = i + 1
#                 amino_acid_seq = translate_sequence(
#                     nucleotide_seq[i:] + Seq("N" * i), translation_table
#                 )
#                 # truncate the stop codon "*" anywhere in the sequence
#                 amino_acid_seq = amino_acid_seq.split("*")[0]
#                 # remove long X sequences at the beginning and end of the sequence
#                 amino_acid_seq = amino_acid_seq.rstrip("X").lstrip("X")
#                 if len(amino_acid_seq) < 5:
#                     amino_acid_seq = Seq("X")
#                 new_header = f"{record.id}_frame_{frame}"
#                 new_record = SeqRecord(amino_acid_seq, id=new_header, description="")
#                 SeqIO.write(new_record, f_out, "fasta")
#             for i in range(3):
#                 frame = i + 4
#                 amino_acid_seq = translate_sequence(
#                     nucleotide_seq.reverse_complement()[i:] + Seq("N" * i),
#                     translation_table,
#                 )
#                 # truncate the stop codon "*" anywhere in the sequence
#                 amino_acid_seq = amino_acid_seq.split("*")[0]
#                 # remove long X sequences at the beginning and end of the sequence
#                 amino_acid_seq = amino_acid_seq.rstrip("X").lstrip("X")
#                 if len(amino_acid_seq) < 5:
#                     amino_acid_seq = Seq("X")
#                 new_header = f"{record.id}_frame_{frame}"
#                 new_record = SeqRecord(amino_acid_seq, id=new_header, description="")
#                 SeqIO.write(new_record, f_out, "fasta")


def translate_fasta(input_file, output_file, translation_table):
    """
    Translates nucleotide sequences in a FASTA file to amino acid sequences
    and writes all 6 reading frames to the output file.
    """
    with open(output_file, "w") as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            nucleotide_seq = Seq(str(record.seq))
            out_seq = Seq("X")
            for i in range(3):
                # frame = i + 1
                amino_acid_seq = translate_sequence(
                    nucleotide_seq[i:] + Seq("N" * i), translation_table
                )
                # truncate the stop codon "*" anywhere in the sequence
                amino_acid_seq = amino_acid_seq.split("*")[0]
                # remove long X sequences at the beginning and end of the sequence
                amino_acid_seq = amino_acid_seq.rstrip("X").lstrip("X")
                if len(amino_acid_seq) < 5:
                    amino_acid_seq = Seq("X")
                if len(amino_acid_seq) > len(out_seq):
                    out_seq = amino_acid_seq
            for i in range(3):
                # frame = i + 4
                amino_acid_seq = translate_sequence(
                    nucleotide_seq.reverse_complement()[i:] + Seq("N" * i),
                    translation_table,
                )
                # truncate the stop codon "*" anywhere in the sequence
                amino_acid_seq = amino_acid_seq.split("*")[0]
                # remove long X sequences at the beginning and end of the sequence
                amino_acid_seq = amino_acid_seq.rstrip("X").lstrip("X")
                if len(amino_acid_seq) < 5:
                    amino_acid_seq = Seq("X")
                if len(amino_acid_seq) > len(out_seq):
                    out_seq = amino_acid_seq
            new_header = f"{record.id}"
            new_record = SeqRecord(out_seq, id=new_header, description="")
            SeqIO.write(new_record, f_out, "fasta")


# Main
def main():
    args = parse_args()
    translate_fasta(args.input_file, args.output_file, args.translation_table)


if __name__ == "__main__":
    main()
