import os
import argparse
import sys
import time
from Bio import Entrez, SeqIO 

# /bin/python /home/alex0414/.vscode-server/extensions/ms-python.python-2024.0.1/pythonFiles/printEnvVariablesToFile.py /home/alex0414/.vscode-server/extensions/ms-python.python-2024.0.1/pythonFiles/deactivate/bash/envVars.txt


# Set your email address for Entrez
def init_Entrez(email):
    Entrez.email = email

# Use Entrez to search for the specified ids and retrieve the list of IDs
def get_batch(refseq_ids):
    handle = Entrez.efetch(db="nucleotide", id=refseq_ids, rettype="fasta", retmode="text")
    sequences = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return sequences

# Create a fasta file for the sequence
def save_to_file(sequence, output_dir, genus):

    # If the output directory doesn't exist, make it
    directory = output_dir + "/" + genus
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Save the genome
    filename = os.path.join(directory, f"{sequence.id.replace(".", "_")}.fa")
    with open(filename, "w") as output_file:
        SeqIO.write(sequence, output_file, "fasta")


# Main
if __name__ == "__main__":

    # Set parameters
    batch_size = 1
    output_dir = "microbial_genomes/filtered_fasta"
    input_file = "NCBI_sequenceIDs.txt"
    genus_filter = ["Streptococcus",
                    "Mycobacterium",
                    "Staphylococcus"]
    
    # Require an email address argument
    parser = argparse.ArgumentParser(description="Download genomes from Refseq, filtered by genus.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    init_Entrez(args.email)
    
    # Keep track of how many of each genome we find
    counts = {genus: 0 for genus in genus_filter}

    # Read in all ids
    with open(input_file, 'r') as file:
        refseq_ids = file.readlines()
    for i in range(len(refseq_ids)):
        refseq_ids[i] = refseq_ids[i].strip()

    # Filter out everything except NC_ ids
    # nc_ids = []
    # for id in refseq_ids:
    #     if "NC_" in id:
    #         nc_ids.append(id)
    # refseq_ids = nc_ids

    # Separate the list into batches
    num_ids = len(refseq_ids)

    # 94200
    for i in range(94197, num_ids, batch_size):
        time.sleep(1)
        sys.stdout.write(f"\rProcessing sequence {i+1}/{num_ids}")
        id_batch = refseq_ids[i:i+batch_size]
        print(id_batch)

        # Get the batch of sequences
        sequences = get_batch(id_batch)
        for sequence in sequences:
            genus = sequence.description.split()[1]

            # If it's a genus we want, save it and add 1 to its count
            if genus in genus_filter:
                save_to_file(sequence, output_dir, genus)
                counts[genus] += 1
    
    # Print results summary
    sys.stdout.write(f"\rProcessed {num_ids} ids, saved:\n")
    for genus in counts:
        print(f"{counts[genus]} {genus} sequences")
