import os
import sys

# This program just takes the fasta files in the input,
# grabs their names, converts their names into an NCBI
# accession version, and then outputs them to a file.

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 2:
    print("Usage: python step_4_2_get_accessions_from_filenames.py <input_file>")
    exit(1)

# Get input and output file path from command-line arguments
input_file = sys.argv[1]

# Check if input_file and output_file are provided
if not input_file:
    print("The input file argument is missing.")
    exit(1)

# Read in the filenames
with open(input_file, 'r') as file:
    file_list = file.readlines()
    file_list = [line.strip() for line in file_list]

# Strip the file extensions
stripped_file_names = [os.path.splitext(file)[0] for file in file_list]

# Replace the last underscore with a period (turn it into an accession version)
# e.g. NC_000001_11 -> NC_000001.11
sequence_ids = []
for file_name in stripped_file_names:
    last_underscore_index = file_name.rfind('_')
    if last_underscore_index != -1:
        modified_string = file_name[:last_underscore_index] + '.' + file_name[last_underscore_index + 1:]
        sequence_ids.append(modified_string)
    else:
        sequence_ids.append(file_name)

# Write the sequence ids to a file
with open(input_file, 'w') as f:
    for i in range(len(sequence_ids)):
        f.write(f"{file_list[i]}\t{sequence_ids[i]}\n")
