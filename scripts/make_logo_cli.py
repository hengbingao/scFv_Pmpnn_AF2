#!/usr/bin/env python

import argparse
from weblogo import *
import os

def parse_fasta(fasta_file):
    """
    Parses a FASTA file and returns a list of sequences.

    Args:
    fasta_file (str): Path to the FASTA file.

    Returns:
    list: A list of sequences (strings).
    """
    sequences = []
    with open(fasta_file, 'r') as f:
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return sequences

def create_sequence_logo(sequences, output_file):
    """
    Creates a sequence logo from a list of sequences and saves it as an image.

    Args:
    sequences (list): A list of protein sequences.
    output_file (str): Path where the sequence logo image will be saved.
    """
    # Convert sequences to a WebLogo SeqList
    seq_list = SeqList(sequences, alphabet=unambiguous_protein_alphabet)
    
    # Create a sequence logo
    logo = LogoData.from_seqs(seq_list)
    logo_options = LogoOptions()
    logo_options.title = "Protein Sequence Logo"
    logo_format = LogoFormat(logo, logo_options)
    
    # Save the sequence logo to a file
    with open(output_file, 'wb') as f:
        f.write(png_formatter(logo, logo_format))

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Generate a sequence logo from a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing protein sequences.")
    parser.add_argument("-o", "--output", help="Path to the output sequence logo image file.", default="sequence_logo.png")

    # Parse command-line arguments
    args = parser.parse_args()

    # Parse the FASTA file to get sequences
    sequences = parse_fasta(args.fasta_file)

    # Create and save the sequence logo
    create_sequence_logo(sequences, args.output)

   # print(f"Sequence logo created and saved to {args.output}")

if __name__ == "__main__":
    main()

