#!/usr/bin/env python3

import argparse
import sys


def convert_fasta_to_phylip(input_file, output_file, name_length=8):
    """
    Convert FASTA format to PAML-compatible PHYLIP format.

    Args:
        input_file: Path to input FASTA file
        output_file: Path to output PHYLIP file
        name_length: Number of characters to keep from sequence names
    """
    try:
        # Read the FASTA file and create PAML-compatible PHYLIP
        with open(input_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    sequences = []
    names = []
    current_seq = ""

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
            # Take first name_length chars after >
            names.append(line[1 : name_length + 1])
            current_seq = ""
        else:
            current_seq += line

    if current_seq:
        sequences.append(current_seq)

    if not sequences:
        print("Error: No sequences found in input file.", file=sys.stderr)
        sys.exit(1)

    try:
        # Write PHYLIP format
        with open(output_file, "w") as f:
            f.write(f" {len(sequences)} {len(sequences[0])}\n")
            for name, seq in zip(names, sequences):
                # Ensure exactly 10 chars (pad with spaces) + 2 spaces + sequence
                padded_name = name.ljust(10)
                f.write(f"{padded_name}  {seq}\n")
    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA format to PAML-compatible PHYLIP format"
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")
    parser.add_argument(
        "-n",
        "--name-length",
        type=int,
        default=8,
        help="Number of characters to keep from sequence names (default: 8)",
    )

    args = parser.parse_args()

    convert_fasta_to_phylip(args.input, args.output, args.name_length)
    print(f"Successfully converted '{args.input}' to '{args.output}'")


if __name__ == "__main__":
    main()
