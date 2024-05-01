"""
translate_satc.py
Daniel Cotter - 05/01/2024

This script takes in a file in SATC (sample, anchor, target, count)
format and outputs a 5 column file with the sample, anchor (nucleotide),
anchor (amino acid), target (amino acid), and count.

# takes input, output file names, and the translation table.
"""

# import statements
import argparse
from Bio.Seq import Seq


# function definitions
def parse_args():
    parser = argparse.ArgumentParser(
        description="This script takes in a file in SATC (sample, anchor, target, count)"
        "format and outputs a 5 column file with the sample, anchor (nucleotide), anchor (amino acid), target (amino acid), and count."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="The input file containing the sample, anchor, target, count",
    )
    parser.add_argument("output_file", type=str, help="The output file")
    parser.add_argument(
        "-t",
        "--translation_table",
        type=str,
        default="Standard",
        help="The translation table to use for translating nucleotide sequences to amino acid sequences",
    )
    parser.add_argument(
        "-r",
        "--try_all_reading_frames",
        action="store_true",
        help="Try all reading frames if stop codons are present in the anchor",
    )
    return parser.parse_args()


def translate_sequence(nucleotide_seq, translation_table):
    """
    Translates a nucleotide sequence to an amino acid sequence
    using the specified translation table
    """
    seq = Seq(nucleotide_seq)
    return seq.translate(table=translation_table)


def main():
    args = parse_args()
    with open(args.input_file, "r") as f_in, open(args.output_file, "w") as f_out:
        for line in f_in:
            sample, anchor, target, count = line.strip().split("\t")
            anchor_nucleotide = anchor
            anchor_amino_acid = translate_sequence(anchor, args.translation_table)
            # if there are stop codons in the anchor, try a different reading frame
            # pad the end of the anchor with Ns to try a different reading frame
            if args.try_all_reading_frames:
                if "*" in anchor_amino_acid:
                    anchor_amino_acid = translate_sequence(
                        anchor[1:] + "N", args.translation_table
                    )
                if "*" in anchor_amino_acid:
                    anchor_amino_acid = translate_sequence(
                        anchor[2:] + "NN", args.translation_table
                    )
            target_amino_acid = translate_sequence(target, args.translation_table)
            if args.try_all_reading_frames:
                if "*" in target_amino_acid:
                    target_amino_acid = translate_sequence(
                        target[1:] + "N", args.translation_table
                    )
                if "*" in target_amino_acid:
                    target_amino_acid = translate_sequence(
                        target[2:] + "NN", args.translation_table
                    )
            f_out.write(
                f"{sample}\t{anchor_nucleotide}\t{anchor_amino_acid}\t{target_amino_acid}\t{count}\n"
            )


if __name__ == "__main__":
    main()
