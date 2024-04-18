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
    # If they included the -f, tell the
    # function to download the full database.
    ./step_4_0_installation.sh -f
else
    # Otherwise, just run the regular, quick install.
    ./step_4_0_installation.sh
fi

# Step 4.1
# Create a list of filenames from the input directory.
./step_4_1_create_filename_list.sh $INPUT_DIRECTORY filenames_and_accessions.tsv

# Step 4.2
# Turn those filenames into accession versions and
# append them to the original file in a new column.
python step_4_2_get_accessions_from_filenames.py filenames_and_accessions.tsv

# Step 4.3
# Use those accession versions to automatically get
# taxids for all that are possible.
# NOTE: Some will fail for various reasons,
# e.g. They were removed, they aren't on NCBI, etc.
./step_4_3_get_taxids_from_accessions.sh filenames_and_accessions.tsv mapping_file.tsv

# Step 4.4
# This takes the original file of filenames and accessions
# and leaves only those that failed to get a taxid.
# If you would like to try again, just repeat steps 3 and 4
# ad nauseum until you feel you won't be able to do better.
python step_4_4_find_missing_accessions.py filenames_and_accessions.tsv mapping_file.tsv

# Step 4.5
# Run the contamination screen. This will produce
# two output files: summary.txt and summary.rpt.
./step_4_5_screen_genomes.sh $INPUT_DIRECTORY mapping_file.tsv $EMAIL_ADDRESS
