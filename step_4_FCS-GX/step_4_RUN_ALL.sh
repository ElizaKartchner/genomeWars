#!/bin/bash

# Check usage
if [ "$#" -ne 2 ] && [ "$#" -ne 3 ]; then
    echo "Usage: $0 <in_dir> <email> [-f]"
    echo "  - in_dir: The input directory containing the fasta files"
    echo "  - email: The user's email, used for the slurm script"
    echo "  -f: Add this flag if you want to download the whole FCS-GX database (470 GB)"
    exit 1
fi

export INPUT_DIRECTORY=$1
export EMAIL_ADDRESS=$2

# Step 4.0: Install packages
if [ "$#" -eq 3 ] && [ "$3" == "-f" ]; then
    ./step_4_0_installation.sh -f
else
    ./step_4_0_installation.sh
fi

# Step 4.1
./step_4_1_create_filename_list.sh $INPUT_DIRECTORY filenames_and_accessions.tsv

# Step 4.2
python step_4_2_get_accessions_from_filenames.py filenames_and_accessions.tsv

# Step 4.3
./step_4_3_get_taxids_from_accessions.sh filenames_and_accessions.tsv mapping_file.tsv

# Step 4.4 (This is only one pass, it won't loop to find any missing taxids)
python step_4_4_find_missing_accessions.py filenames_and_accessions.tsv mapping_file.tsv

# Step 4.5
./step_4_5_screen_genomes.sh $INPUT_DIRECTORY mapping_file.tsv $EMAIL_ADDRESS
