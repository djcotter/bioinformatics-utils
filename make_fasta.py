"""
make_fasta.py
Daniel Cotter -- 04/08/24

Create a fasta file from a list of sequences provided from standard input
"""

import sys


# Main function
def main():
    sequences = []
    for line in sys.stdin:
        if line.startswith(">"):
            continue
        sequences.append(line.strip())
    # write the seqs to standard output
    for i, seq in enumerate(sequences):
        print(f">seq_{i}")
        print(seq)
    return None


# run main function
if __name__ == "__main__":
    main()
