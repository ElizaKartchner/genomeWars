import argparse

# This program checks to see which genomes got missed
# when grabbing taxids from the database en masse.

# Create the parser
parser = argparse.ArgumentParser(description="Process some files.")

# Add the arguments
parser.add_argument('all_accessions', type=str, help='The all_accessions file')
parser.add_argument('taxid_map', type=str, help='The taxid_map file')
parser.add_argument('-o', '--output', type=str, help='The output file for the missing accessions')

# Parse the arguments
args = parser.parse_args()

# Get filepaths from command-line arguments
all_accessions_file = args.all_accessions
mapping_file = args.taxid_map

# If the output file argument is not provided, default to all_accessions_file
output_file_path = args.output if args.output else all_accessions_file

# Check if files are provided
if not all_accessions_file:
    print("The accession list argument is missing.")
    exit(1)
if not mapping_file:
    print("The mapping file argument is missing.")
    exit(1)
if not output_file_path:
    print("The output file argument is missing.")
    exit(1)

# Read in the files
with open(all_accessions_file, 'r') as file:
    all_accessions = file.readlines()
    all_accessions = [line.strip() for line in all_accessions]
    all_accessions = [line.split() for line in all_accessions]
with open(mapping_file, 'r') as file:
    taxid_map = file.readlines()
    taxid_map = [line.strip() for line in taxid_map]
    taxid_map = [line.split() for line in taxid_map]

# Get the list of accession versions that got
# lost in the process of getting the taxids
missed_accessions = []
for line in all_accessions:
    found = False
    for map_line in taxid_map:
        if line[0] == map_line[0]:
            found = True
            break
    if not found:
        missed_accessions.append(line)

# Write the lost accessions to a file
with open(output_file_path, 'w') as f:
    for accession in missed_accessions:
        f.write(r"{'\t'.join(accession)}")
        f.write('\n')
        