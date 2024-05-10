#!/usr/bin/env python

import argparse
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
import os

#input_dir = 'input_dir/HA_bound/'

def plddt_ptm_plots(input_dir, prefix=""):
    home = os.getcwd()
    os.chdir(input_dir)
    json_paths = glob.glob('structures/*rank_001*.json')
    # structures/HA_bound_seq_WT_scores_rank_001_alphafold2_ptm_model_4_seed_000.json
    
    all_plddt = []
    all_ptm   = []
    all_names = []
    for jpath in json_paths:
        with open(jpath) as file:
            data = json.load(file)
    
            plddt = data['plddt']
            ptm   = data['ptm']
            all_plddt.append( round(np.mean(plddt), 2) )
            all_ptm.append( ptm )
            idx = jpath.find('seq')
            if jpath[idx + 3] == '_':
                all_names.append('WT')
            
            else:
                #  HA_bound_seq14_score0.
                start_string = jpath[idx:idx+6]
                ss_list      = start_string.split('_')
                num          = ss_list[0][3:]
                all_names.append(num)
    
    #print(all_names)
    #print(all_plddt)
    #print(all_ptm)
    
    combo           = zip(all_plddt, all_ptm, all_names)
    sorted_combined = sorted(combo, key=lambda x: x[0])
    primary_sorted, second_sorted, third_sorted = zip(*sorted_combined)
    plddt_sorted    = list(primary_sorted)
    ptm_sorted      = list(second_sorted)
    names_sorted    = list(third_sorted)
    ymin            = min(plddt_sorted) - 5
    ymax            = max(plddt_sorted) + 5

    plt.figure(1)
    fig, ax = plt.subplots(1)
    bars = ax.bar(names_sorted, plddt_sorted, color='skyblue')
    
    highlight = 'WT'
    labels = names_sorted
    for bar, label in zip(bars, labels):
        if label == highlight:
            bar.set_color('#FF0000')
    
    
    # Adding title and labels
    ax.set_title('pLDDT vs ProteinMPNN seq')
    ax.set_xlabel('seq ID')
    ax.set_ylabel('pLDDT')
    ax.set_ylim(ymin, ymax)
    fig.savefig(f'{prefix}plddt.png', dpi=600)
    
    
    ####################
    combo           = zip(all_ptm, all_plddt, all_names)
    sorted_combined = sorted(combo, key=lambda x: x[0])
    primary_sorted, second_sorted, third_sorted = zip(*sorted_combined)
    plddt_sorted    = list(second_sorted)
    ptm_sorted      = list(primary_sorted)
    names_sorted    = list(third_sorted)
    ymin            = min(ptm_sorted) - .05
    ymax            = max(ptm_sorted) + .05


    
    
    plt.plot(2)
    fig, ax = plt.subplots(1)
    bars = ax.bar(names_sorted, ptm_sorted, color='skyblue')
    
    labels = names_sorted
    highlight = 'WT'
    for bar, label in zip(bars, labels):
        if label == highlight:
            bar.set_color('#FF0000')
    
    
    # Adding title and labels
    ax.set_title('pTM vs ProteinMPNN seq')
    ax.set_xlabel('seq ID')
    ax.set_ylabel('pTM')
    ax.set_ylim(ymin, ymax)
    fig.savefig(f'{prefix}ptm.png', dpi=600)
    os.chdir(home)

def main():
   # Set up the argument parser
    parser = argparse.ArgumentParser(description="Point to a dir with colabfold_batch structures to generate plddt/ptm plots")
    parser.add_argument("dir", help="Input directory to be analyzed")
    parser.add_argument("--prefix", default="", required=False, help="An optional prefix for the output ptm/plddt plots for organization")

    # Parse command-line arguments
    args = parser.parse_args()
    
    if args.prefix != "" and args.prefix[-1] != '_':
        args.prefix += '_'
    plddt_ptm_plots(args.dir, args.prefix)


if __name__ == "__main__":
    main()

