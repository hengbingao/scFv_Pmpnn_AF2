#!/usr/bin/env python

import argparse
import os

def simplify_fasta(fa_file, new_seq=''):
    def process_file(file_path, new_seq):
        # Extracting the descriptor from the file name
        descriptor_base = os.path.basename(file_path).split('.fa')[0].strip('_')

        amended_sequences = []
        seq_counter = 0
        for line in open(file_path, 'r'):
            if line.startswith(">"):
                # Extract score if it exists and isn't the first sequence
                score = ""
                if "score=" in line and seq_counter != 0:
                    score_value = line.split("score=")[1].split()[0]  # assuming score= is followed by space or end of line
                    score = f"_score{score_value[:-1]}"
                if seq_counter == 0:
                    #seq_counter += 1
                    seq_counter = '_WT'
                amended_sequences.append(f">{descriptor_base}_seq" + str(seq_counter) + score + "\n")
                if seq_counter == '_WT':
                    seq_counter = 0
                seq_counter += 1
            else:
                if new_seq == '':
                    amended_sequences.append(line)
                else:
                    line = line.strip() + ":" + new_seq + "\n"
                    amended_sequences.append(line)
        return amended_sequences

    # Check if fa_file is a directory
    if os.path.isdir(fa_file):
        all_sequences = []
        for file_name in os.listdir(fa_file):
            if file_name.endswith('.fa') or file_name.endswith('.fasta'):
                file_path = os.path.join(fa_file, file_name)
                all_sequences.extend(process_file(file_path, new_seq))
        return all_sequences
    else:
        return process_file(fa_file, new_seq)

#modified_fa = simplify_fasta('output_dir/seqs/mbg17_bound.fa')
#with open("simplified_output.fasta", "w") as out_file:
#    out_file.writelines(modified_fa)

def main():

    parser = argparse.ArgumentParser(description="Simplify a ProteinMPNN output fasta/fa file")
    parser.add_argument("fa_file", type=str,  help="Path to the .fa/.fasta file.")
    
    args        = parser.parse_args()
    simple_list = args.fa_file.split('.')
    simple      = '.'.join(simple_list[:-1]) + '.fa'
    simpler     = simple.split('/')[-1]

    print(f'writing file: {simple}')
    modified_fa = simplify_fasta(args.fa_file)
    with open(simpler, "w") as out_file:
        out_file.writelines(modified_fa)

if __name__ == "__main__":
    main()

