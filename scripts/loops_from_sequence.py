#!/usr/bin/env python

from anarci import anarci
import sys
import argparse
import numpy as np
from scipy.spatial.distance import cdist
from collections import defaultdict


def loops_from_sequence(input_seq, numbering='kabat'):
    sequences = [("id1", input_seq)] 

    # Run ANARCI
    # this `results` product is the *worst* data structure I have ever seen in my life -- JD
    results = anarci(sequences, scheme=numbering)
    
    # Kabat numbering for Heavy Chain CDRs
    # numbers for loops citation: https://www.novoprolabs.com/tools/cdr
    cdr_h1 = [str(i) for i in range(31, 36)]
    cdr_h2 = [str(i) for i in range(50, 66)]
    cdr_h3 = [str(i) for i in range(95, 103)]
    
    # Kabat numbering for Light Chain CDRs
    cdr_l1 = [str(i) for i in range(24, 35)]
    cdr_l2 = [str(i) for i in range(50, 57)]
    cdr_l3 = [str(i) for i in range(89, 98)]
    
    # Define all CDRs in a dictionary for easy reference
    cdr_loops = {
        'H': {
            'CDR H1': cdr_h1,
            'CDR H2': cdr_h2,
            'CDR H3': cdr_h3
        },
        'L': {
            'CDR L1': cdr_l1,
            'CDR L2': cdr_l2,
            'CDR L3': cdr_l3
        }
    }
    

    all_loops  = []
    num_chains = len(results[0][0])
    
    for chain_number in range(num_chains):
        highlights = results[0][0][chain_number][0]
        chain_type = results[1][0][chain_number]['chain_type']

        if chain_type.upper() == 'K':
            chain_type = 'L'
        elif chain_type.upper() in ['A', 'D', 'E', 'G', 'M']:
            chain_type = 'H'
    
        # data structure for storing extracted/parsed info 
        loop_seqs = {
                'L1' : ['', []],
                'L2' : ['', []],
                'L3' : ['', []]
        }
        
        # Process each highlighted residue
        for idx, hl in enumerate(highlights):
    
            (num, ins), res = hl
            full_num        = str(num)

            # Check if the residue is part of any CDR loop for the specified chain type
            cdr_found = False
            for cdr_name, cdr_range in cdr_loops[chain_type].items():
                if full_num in cdr_range:
                    cdr_found = True
                    loop_seqs[f'L{cdr_name[-1]}'][0] += res 
                    loop_seqs[f'L{cdr_name[-1]}'][1].append(idx) 
                    break
                
        
        for loop, data in loop_seqs.items():
            seq, resis = data
            seq_n = len(seq)
            start = input_seq.find(seq)

            for i in range(start, start+seq_n):
                all_loops.append(i)
    

    return all_loops


def pdb_to_sequence_manual(pdb_path):
    """Manually extracts the amino acid sequence from a PDB file."""
    # Define a dictionary to convert three-letter codes to one-letter codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U', 'PYL': 'O'  # Including selenocysteine (Sec) and pyrrolysine (Pyl)
    }

    sequences = {}
    last_residue_id = None

    try:
        with open(pdb_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    residue_name = line[17:20].strip()
                    residue_id = line[22:27].strip()  # Including insertion code

                    if chain_id not in sequences:
                        sequences[chain_id] = ''

                    if residue_id != last_residue_id:
                        if residue_name in three_to_one:
                            sequences[chain_id] += three_to_one[residue_name]
                        else:
                            sequences[chain_id] += '?'  # Unrecognized residue
                    last_residue_id = residue_id

    except Exception as e:
        return f"Error reading PDB file: {e}"

    return sequences


def parse_pdb(pdb_path):
    """ Parse PDB to extract coordinates and other data per atom. """
    atoms = []
    with open(pdb_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_data = {
                    'chain': line[21],
                    'res_id': line[22:26].strip() + line[21],  # Combining residue number with chain identifier directly
                    'res_name': line[17:20].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54])
                }
                atoms.append(atom_data)
    return atoms


# Function to parse the PDB file and extract chain information
def parse_pdb_file(file_path):
    chains = defaultdict(list)
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                chain_id = line[21]
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                chains[chain_id].append((residue_number, (x, y, z)))
    return chains


def fill_gaps_and_remove_isolated_residues(residues):
    # Start with the original list of residues in contact
    filled_residues = sorted(set(residues))

    # Initialize a list to hold the intermediate set of residues, including gap-filled ones
    intermediate_residues = []

    # Go through the sorted list and fill in the gaps
    i = 0
    while i < len(filled_residues) - 1:
        intermediate_residues.append(filled_residues[i])

        # Check the gap between the current and next residue
        gap = filled_residues[i + 1] - filled_residues[i]

        # If the gap is 2 or 3 (indicating 1 or 2 residues missing between them), fill it
        if gap <= 3:
            # Add the missing residues to the intermediate list
            intermediate_residues.extend(range(filled_residues[i] + 1, filled_residues[i + 1]))

        i += 1

    # Add the last residue since it's not covered in the loop
    intermediate_residues.append(filled_residues[-1])

    # Remove isolated residues
    final_residues = []
    for res in intermediate_residues:
        # Check if the residue has neighbors within 5 residue numbers
        has_neighbors = any(abs(res - other) <= 3 for other in intermediate_residues if other != res)
        if has_neighbors:
            final_residues.append(res)

    return sorted(final_residues)


def find_close_residues_in_longest_chain(file_path, input_residue_ids, dist=4):
    chains = parse_pdb_file(file_path)
    # Identify the longest chain
    longest_chain_id = max(chains.keys(), key=lambda k: len(set([res[0] for res in chains[k]])))
    longest_chain = chains[longest_chain_id]

    # Collect positions of input residues if they are in the longest chain
    input_positions = [atom_pos for res_num, atom_pos in longest_chain if res_num in input_residue_ids]

    # Find residues in the longest chain close to any atom in the input list
    close_residues = set()
    for atom_pos in input_positions:
        for res_num, atom_pos_l in longest_chain:
            if calculate_distance(atom_pos, atom_pos_l) <= dist:
                close_residues.add(res_num)

    # Optionally process the results to fill gaps and remove isolated residues, if required
    return fill_gaps_and_remove_isolated_residues(list(close_residues)), longest_chain_id


'''
def find_close_residues(pdb_path, input_residues):
    atoms = parse_pdb(pdb_path)

    # Organize atoms by chains and count chain lengths
    chains = {}
    for atom in atoms:
        chain_key = atom['chain']
        if chain_key not in chains:
            chains[chain_key] = []
        chains[chain_key].append(atom)

    # Determine the longest chain
    longest_chain = max(chains, key=lambda k: len(chains[k]))

    # Collect all atom coordinates and residues in the longest chain
    chain_atoms = chains[longest_chain]
    chain_residues = {atom['res_id']: (atom['res_name'], np.array([atom['x'], atom['y'], atom['z']])) for atom in chain_atoms}

    input_coords = []
    chain_coords = []
    chain_res_keys = []

    for res_id, res_data in chain_residues.items():
        chain_coords.append(res_data[1])
        chain_res_keys.append(res_id)
        if res_id in input_residues:
            input_coords.append(res_data[1])

    if not input_coords:
        return "No input residues found in the longest chain."

    distances = distance.cdist(input_coords, chain_coords)
    close_residue_ids = set()

    for i, row in enumerate(distances):
        close_residues = [chain_res_keys[j] for j, d in enumerate(row) if d < 4]
        close_residue_ids.update(close_residues)

    # Extract numbers and return a sorted list of integers
    close_residue_numbers = [int(res_id[:-1]) for res_id in close_residue_ids if res_id.endswith('A')]
    input_residue_numbers = [int(res[:-1]) for res in input_residues if res.endswith('A')]
    combined_residue_numbers = list(set(close_residue_numbers + input_residue_numbers))
    combined_residue_numbers.sort()

    return combined_residue_numbers
'''


def find_close_residues(pdb_file_path, residues, cutoff=4.0):
    """
    Find neighboring residues within a given cutoff distance from a list of input residues.

    Args:
        pdb_file_path (str): Path to the PDB file.
        residues (list): List of residue IDs (e.g., ['47A']) to find neighbors for.
        cutoff (float): Maximum distance (in Angstroms) to consider a residue as a neighbor.

    Returns:
        list: List of unique residue numbers that are within the cutoff distance of any residue in the input list.
    """
    coords = []
    res_ids = []
    chains = []

    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                chain = line[21]
                res_id = int(line[22:26].strip())
                res_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coords.append([x, y, z])
                res_ids.append(res_id)
                chains.append(chain)

    coords = np.array(coords)
    res_ids = np.array(res_ids)
    chains = np.array(chains)

    input_chains = [res[-1] for res in residues]
    input_res_ids = [int(res[:-1]) for res in residues]

    neighboring_residues = set()

    for res_id, chain in zip(input_res_ids, input_chains):
        res_mask = (res_ids == res_id) & (chains == chain)
        res_coords = coords[res_mask]
        if res_coords.ndim == 1:
            res_coords = res_coords.reshape(1, -1)

        other_res_mask = ~res_mask
        other_res_coords = coords[other_res_mask]
        other_res_ids = res_ids[other_res_mask]
        other_chains = chains[other_res_mask]

        if other_res_coords.ndim == 1:
            other_res_coords = other_res_coords.reshape(1, -1)

        pairwise_dists = cdist(res_coords, other_res_coords)
        neighboring_residues.update([other_res_id for other_res_id, other_chain, dists in zip(other_res_ids, other_chains, pairwise_dists.T) if np.any(dists <= cutoff) and other_chain in input_chains])

    output_residues = [res_id for res_id in neighboring_residues]

    return sorted(output_residues)


def get_longest_chain_length(pdb_path):
    """Returns the length of the longest chain in a PDB file by counting unique residues."""
    chain_residues = {}

    try:
        with open(pdb_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):  # Ensures we only process ATOM lines
                    chain_id = line[21]
                    residue_id = line[22:26].strip() + line[26]  # Combining residue number with insertion code

                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = set()
                    
                    chain_residues[chain_id].add(residue_id)
    
    except FileNotFoundError:
        return "PDB file not found."

    # Find the longest chain by comparing the size of residue sets
    longest_chain_length = max((len(residues) for residues in chain_residues.values()), default=0)

    return longest_chain_length


def main():

    # I haven't tested it, but I almost promise this breaks if input is anything other than a PDB
    parser = argparse.ArgumentParser(description="returns an index of the CDR loops as identified by ANARCI.")
    parser.add_argument("fasta_file", type=str,                  help="Path to the fasta file to identify CDR loops")
    parser.add_argument("--verbose",  action='store_true',       help="True/False: print PyMOL selection logic.")
    parser.add_argument("--scheme",   type=str, default="kabat", help="How do we want to number the antibody H/L chains? Default = kabat.")
    parser.add_argument("--output", choices=['loops', 'lss', 'framework'],
                    help="Output option: loops, lss (loops AND secondary shell [aka Vernier zone region]), or framework.", default='framework')


    args = parser.parse_args()
    
    seqs = []
    if args.fasta_file.split('.')[-1] == 'fasta' or args.fasta_file.split('.')[-1] == 'fa':
        with open(args.fasta_file) as file:
            for line in file:
                if line[0] != '>' and len(line) > 5: # find only protein sequences in fasta file
                    seqs.append(line.replace('\n', ''))
    
    elif args.fasta_file.split('.')[-1] == 'pdb':
        sequences = pdb_to_sequence_manual(args.fasta_file)

        for chain, seq in sequences.items():
            if len(seq) > 30: # try and only extract VH/VL chains and NOT any peptides
                seqs.append(seq)

    all_loops = []
    # sequence-wise find all loops
    for seq in seqs:
        all_loops += loops_from_sequence(seq, numbering=args.scheme)

    # sequence-wise find secondary shell (unformally ID the vernier zone region)
    if args.fasta_file.split('.')[-1] == 'pdb':

        a = []
        # extra here
        [a.append(f'{x}A') for x in all_loops]    
        residues = find_close_residues(args.fasta_file, a, cutoff=4)
        vernier_and_loops = fill_gaps_and_remove_isolated_residues(residues)
        length   = get_longest_chain_length(args.fasta_file)
        all_resi = list(range(1, length+1))

    framework = [x for x in all_resi if x not in vernier_and_loops]

    if   args.output == 'loops':
        out = all_loops
    elif args.output == 'lss':
        out = vernier_and_loops
    elif args.output == 'framework':
        out = framework

    if args.verbose:
        a = "+".join(str(x) for x in all_loops)
        print('### some potentially helpful PyMOL selection logic ###')
        print('Loops:')
        print(f'select resi {a}\n')

        a = "+".join(str(x) for x in vernier_and_loops)
        print('Loops and secondary shell:')
        print(f'select resi {a}\n')

        a = "+".join(str(x) for x in framework)
        print('framework:')
        print(f'select resi {a}\n')
        print('######################################################')
    a = " ".join(str(x) for x in out)
    print(a)
    return out
if __name__ == "__main__":
    main()
