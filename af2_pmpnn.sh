#!/bin/bash
# @jderoo modified for parallel processing — improved version

# ================================
# Set ColabFold Cache Directory
# ================================
export COLABFOLD_CACHE_DIR="/group/ll010/hgao/apps/colabfold_cache"
export XDG_CACHE_HOME="/group/ll010/hgao/apps/colabfold_cache"

##################
### USER INPUT ###
##################
if [ -z "$1" ]; then
    echo "Usage: bash af2_pmpnn_parella.sh <input_folder>"
    exit 1
fi

folder_with_pdbs="$1"
seqs_per_run=15
determine_CDRs='martin'
ss_near_CDRs=3
simple_grab=true
output_dir="$folder_with_pdbs/out"
to_design="framework"
linker_seq=""
PMPNN='/group/ll010/hgao/apps/ProteinMPNN'

##################
### SETUP PATHS ###
##################
mkdir -p "$output_dir"
mkdir -p "$output_dir/seqs"

pdb_count=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.pdb" | wc -l)
fasta_count=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.fasta" | wc -l)
pdbs=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.pdb")
fastas=$(find "$folder_with_pdbs" -maxdepth 1 -type f -name "*.fasta")

echo ">>> Processing folder: $folder_with_pdbs"
echo ">>> Found $pdb_count PDBs and $fasta_count FASTAs"

# ===========================================
# Step 1: Move PDBs into subfolders
# ===========================================
if [ "$pdb_count" -gt 0 ]; then
    for pdb_file in $pdbs; do
        base_name=$(basename "$pdb_file" .pdb)
        if [[ "$base_name" == *"seq"* ]]; then
            echo "Error: 'seq' pattern found in PDB filename: $base_name"
            exit 1
        fi
        mkdir -p "$folder_with_pdbs/$base_name"
        mv "$pdb_file" "$folder_with_pdbs/$base_name"
    done
fi
echo ">>> Finished moving PDB files."

# ===========================================
# Step 2: Process FASTAs → structures
# ===========================================
if [ "$fasta_count" -gt 0 ]; then
    for fasta_file in $fastas; do
        base_name=$(basename "$fasta_file" .fasta)
        mkdir -p "$folder_with_pdbs/$base_name"

        if [[ "$base_name" == *"_scfv" ]]; then
            echo ">>> [$base_name] Detected scFv format — skipping scfv_anarci.py conversion."
            fasta_to_use="$fasta_file"
        else
            echo ">>> [$base_name] Detected raw VH/VL fasta — converting to scFv..."
            /group/ll010/hgao/apps/scFv_Pmpnn_AF2/scripts/scfv_anarci.py "$fasta_file" --output "${base_name}_scfv.fasta" --polyG
            fasta_to_use="${base_name}_scfv.fasta"
        fi

        echo ">>> [$base_name] Predicting initial structure with AlphaFold2..."
        colabfold_batch --templates --num-recycle 6 --model-type alphafold2_multimer_v2 "$fasta_to_use" "${base_name}_initial_structures"

        best_pdb=$(ls "${base_name}_initial_structures"/*rank_001*pdb 2>/dev/null | head -n1)
        if [ -z "$best_pdb" ]; then
            echo "Error: No PDB generated for $base_name"
            exit 1
        fi

        cp "$best_pdb" "./${base_name}.pdb"
        mv "./${base_name}.pdb" "$folder_with_pdbs/$base_name"
        mv "$fasta_to_use" "$folder_with_pdbs/$base_name"
        mv "$fasta_file" "$folder_with_pdbs/$base_name"
        rm -rf "${base_name}_initial_structures"
    done
fi
echo ">>> Finished generating structures."

# ===========================================
# Step 3: Main loop for each subdirectory
# ===========================================
for dir in "$folder_with_pdbs"/*; do
    if [ ! -d "$dir" ] || [[ "$dir" == "$output_dir" ]]; then
        continue
    fi

    dir=${dir%/}
    echo ">>> Processing directory: $dir"

    pdb_file=$(find "$dir" -type f -name "*.pdb" | head -n1)
    if [ -z "$pdb_file" ]; then
        echo "Warning: No PDB found in $dir, skipping."
        continue
    fi

    trimmed=$(basename "$dir")

    # Output files made unique to prevent parallel conflicts
    path_for_parsed_chains="$output_dir/parsed_pdbs_${trimmed}.jsonl"
    path_for_assigned_chains="$output_dir/assigned_pdbs_${trimmed}.jsonl"
    path_for_fixed_positions="$output_dir/fixed_pdbs_${trimmed}.jsonl"

    chains_to_design=$(/group/ll010/hgao/apps/scFv_Pmpnn_AF2/scripts/longest_chain.py "$pdb_file")

    if [ "$determine_CDRs" = "structure" ]; then
        IFS=" " read -r design_only_positions <<< $(/group/ll010/hgao/apps/scFv_Pmpnn_AF2/scripts/find_loops.py "$pdb_file" --output "$to_design")
    else
        cmd="/group/ll010/hgao/apps/scFv_Pmpnn_AF2/scripts/loops_from_sequence.py $pdb_file --scheme $determine_CDRs --output $to_design --dist $ss_near_CDRs --simple_grab $simple_grab"
        [ -n "$linker_seq" ] && cmd+=" --linker-seq $linker_seq"
        IFS=" " read -r design_only_positions <<< $($cmd)
    fi

    echo ">>> [$trimmed] Running ProteinMPNN design..."
    python "$PMPNN/helper_scripts/parse_multiple_chains.py" --input_path="$dir" --output_path="$path_for_parsed_chains"
    python "$PMPNN/helper_scripts/assign_fixed_chains.py" --input_path="$path_for_parsed_chains" --output_path="$path_for_assigned_chains" --chain_list "$chains_to_design"
    python "$PMPNN/helper_scripts/make_fixed_positions_dict.py" --input_path="$path_for_parsed_chains" --output_path="$path_for_fixed_positions" --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

    python "$PMPNN/protein_mpnn_run.py" \
        --jsonl_path "$path_for_parsed_chains" \
        --chain_id_jsonl "$path_for_assigned_chains" \
        --fixed_positions_jsonl "$path_for_fixed_positions" \
        --out_folder "$output_dir" \
        --num_seq_per_target "$seqs_per_run" \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1 \
        --use_soluble_model

    echo ">>> [$trimmed] Simplifying designed fasta..."
    /group/ll010/hgao/apps/scFv_Pmpnn_AF2/scripts/simplify_fasta.py "$output_dir/seqs/${trimmed}.fa"

    tail -n +3 "$output_dir/seqs/${trimmed}.fa" > "$dir/${trimmed}_noWT.fa"
    head -n +2 "$output_dir/seqs/${trimmed}.fa" > "$dir/${trimmed}_WT.fa"

    echo ">>> [$trimmed] Running AlphaFold2 on redesigned sequence..."
    colabfold_batch --templates --num-recycle 6 --model-type alphafold2_multimer_v2 "$dir/${trimmed}_noWT.fa" "$dir/af2_output"

    echo ">>> [$trimmed] Completed successfully."
done

echo "? All jobs finished for $folder_with_pdbs!"
