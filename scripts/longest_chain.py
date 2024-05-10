#!/usr/bin/env python3

import argparse
import sys

def read_pdb(file_path):
    """
    Reads a PDB file and returns a dictionary with chain IDs as keys and the count of unique residues in each chain as values.
    """
    chains = {}
    residues = set()
    try:
        with open(file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chain_id = line[21]
                    residue_number = line[22:26].strip()
                    residue_key = f"{chain_id}-{residue_number}"
                    if residue_key not in residues:
                        residues.add(residue_key)
                        if chain_id in chains:
                            chains[chain_id] += 1
                        else:
                            chains[chain_id] = 1
    except IOError:
        print(f"Error: Could not read file {file_path}")
        sys.exit(1)
    return chains

def find_longest_chain(chains):
    """
    Finds and returns the chain with the most residues.
    """
    return max(chains, key=chains.get)

def main():
    parser = argparse.ArgumentParser(description='Find the longest chain in a protein structure based on residue count.')
    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
    args = parser.parse_args()

    chains = read_pdb(args.pdb_file)
    
    if chains:
        longest_chain = find_longest_chain(chains)
        print(longest_chain)
    else:
        print("No chains found in the provided PDB file.")

if __name__ == "__main__":
    main()

