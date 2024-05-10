#!/usr/bin/env python
"""
Created on Sat Feb 10 04:20:08 2024

@author: jderoo
"""

from collections import defaultdict
import math
from scipy.spatial.distance import cdist
import numpy as np
import argparse

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




# Function to calculate distance between two atoms
def calculate_distance(atom1, atom2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(atom1, atom2)))



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
        has_neighbors = any(abs(res - other) <= 5 for other in intermediate_residues if other != res)
        if has_neighbors:
            final_residues.append(res)

    return sorted(final_residues)



# Function to find residues in the longest chain within 5 Angstroms of the shortest chain
def find_contacting_residues(file_path, dist=6):
    chains = parse_pdb_file(file_path)
    # Identify the smallest and largest chains
    smallest_chain_id = min(chains.keys(), key=lambda k: len(set([res[0] for res in chains[k]])))
    largest_chain_id  = max(chains.keys(), key=lambda k: len(set([res[0] for res in chains[k]])))
    smallest_chain    = chains[smallest_chain_id]
    largest_chain     = chains[largest_chain_id]

    # Find residues in the largest chain close to any atom in the smallest chain
    close_residues = set()
    for _, atom_pos in smallest_chain:
        for res_num, atom_pos_l in largest_chain:
            if calculate_distance(atom_pos, atom_pos_l) <= dist:
                close_residues.add(res_num)

    return fill_gaps_and_remove_isolated_residues(close_residues), largest_chain_id



def parse_pdb_for_atoms_specific_chain(file_path, chain_id):
    all_atoms = defaultdict(lambda: [])
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[21] == chain_id:
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                all_atoms[residue_number].append([x, y, z])
    # Convert lists to numpy arrays for efficient distance calculations
    for res_id in all_atoms:
        all_atoms[res_id] = np.array(all_atoms[res_id])
    return all_atoms

def find_second_shell(pdb_file, chain_id, input_residue_list, dist=4):
    # Parse the PDB file for atom coordinates in the specified chain
    all_atoms = parse_pdb_for_atoms_specific_chain(pdb_file, chain_id)

    # Extract atom coordinates for input residues
    target_coords = np.vstack([all_atoms[res] for res in input_residue_list if res in all_atoms])

    close_residues = set()
    # Calculate distances from target residues to all other residues in the chain
    for res, coords in all_atoms.items():
        if res in input_residue_list:
            continue  # Skip the residues already in the input list
        distances = cdist(target_coords, coords)
        if np.any(distances <= dist):
            close_residues.add(res)

    return fill_gaps_and_remove_isolated_residues(close_residues)



def find_loops(file_path, verbose=False):

    contacting_residues, scFv_chain  = find_contacting_residues(file_path, dist=6)
    second_shell                     = find_second_shell(file_path, scFv_chain, contacting_residues, dist=4)
    second_shell_only                = [x for x in second_shell if x not in contacting_residues]
    pymol                            = '+'.join([str(x) for x in contacting_residues])
    pymol_sso                        = '+'.join([str(x) for x in second_shell_only])


    scFv_length = parse_pdb_file(file_path)[scFv_chain][-1][0]
    keepers     = second_shell_only + contacting_residues
    mutables    = [int(x) for x in range(1,scFv_length+1) if x not in keepers]
    pymol_rest  = '+'.join([str(x) for x in range(1,scFv_length+1) if x not in keepers])

    if verbose:
        print('\n############### Direct Contact ###################')
        print(f'select chain {scFv_chain} and resi {pymol}')
        print('--------------------------------------------------')

        print('################ Second Shell ####################')
        print(f'select chain {scFv_chain} and resi {pymol_sso}')
        print('--------------------------------------------------')

        print('################  Framework  #####################')
        print(f'select chain {scFv_chain} and resi {pymol_rest}')
        print('--------------------------------------------------')

    return contacting_residues, second_shell_only, mutables

def main():
    parser = argparse.ArgumentParser(description="Find loops and framework of an scFv structure. If this is called on the command line, the INVERSE is returned by default for easy input into ProteinMPNN - the *FRAMEWORK* is returned, NOT the loops! See --output flag for behavior.")
    parser.add_argument("pdb_file", type=str,  help="Path to the PDB file.")
    parser.add_argument("--verbose", action='store_true', help="True/False: print PyMOL selection logic.")
    parser.add_argument("--output", choices=['loops', 'ss', 'lss', 'framework'],
                    help="Output option: loops, ss (secondary shell), lss (loops AND secondary shell), or framework.", default='framework')

    args = parser.parse_args()

    loops, ss, framework = find_loops(args.pdb_file, verbose=args.verbose)
    combo = loops + ss

    if args.output   == 'loops':
        out   = ' '.join([str(x) for x in loops])
    elif args.output == 'ss':
        out   = ' '.join([str(x) for x in ss])
    elif args.output == 'lss':
        out   = ' '.join([str(x) for x in combo])
    elif args.output == 'framework':
        out   = ' '.join([str(x) for x in framework])
    else:
        print('flag error! please use an appropriate flag (loops, ss, lss, framework) for --output')

    print(out)




if __name__ == "__main__":
    main()
