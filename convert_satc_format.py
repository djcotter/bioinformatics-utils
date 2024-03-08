"""
convert_satc_format.py
Daniel Cotter -- 02/16/24

Convert satc dump files from the old format (count, anchor+target, sample)
to the new format (sample, anchor, target, count)
"""

import argparse
import csv


# Parse commandline arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert satc dump files to new format"
    )
    parser.add_argument("input", help="Input file", type=str)
    parser.add_argument("output", help="Output file", type=str)
    parser.add_argument("-k", "--anchor", help="Anchor length", default=27, type=int)
    return parser.parse_args()


# Process each line of the input file and write to the output file
def process_file(input_file: str, output_file: str, k: int) -> None:
    with open(input_file, "r") as f:
        with open(output_file, "w") as o:
            writer = csv.writer(o, delimiter="\t")
            for line in f:
                line = line.strip().split()
                count = line[0]
                anchor = line[1][:k]
                target = line[1][k:]
                sample = line[2]
                writer.writerow([sample, anchor, target, count])
    return None


# Main function
def main():
    args = parse_args()
    process_file(args.input, args.output, args.anchor)


# run main function
if __name__ == "__main__":
    main()
