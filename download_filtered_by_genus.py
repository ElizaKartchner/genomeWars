import os
import argparse
import sys
import time
import pdb
from Bio import Entrez, SeqIO 

# /bin/python /home/alex0414/.vscode-server/extensions/ms-python.python-2024.0.1/pythonFiles/printEnvVariablesToFile.py /home/alex0414/.vscode-server/extensions/ms-python.python-2024.0.1/pythonFiles/deactivate/bash/envVars.txt

# Set your email address for Entrez
def init_Entrez(email):
    Entrez.email = email

# Use Entrez to search by the accession number and retrieve the list of fasta files
def get_file_batch(accession_nums):
    handle = Entrez.efetch(db="nucleotide", id=accession_nums, rettype="fasta", retmode="text")
    sequences = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return sequences

# # Use Entrez to search by the accession number and retrieve summaries which contain the taxid
# def get_summary_batch(accession_nums):
#     # Turn the list into a single, comma-separated string of values
#     accession_nums = ",".join(accession_nums)
#     handle = Entrez.esummary(db="nuccore", id=accession_nums, retmode="xml")
#     summaries = Entrez.read(handle)
#     handle.close()
#     return summaries

# # Create a fasta file for the sequence
# def save_to_file(sequence, output_dir, genus, tax_id):

#     # If the output directory doesn't exist, make it
#     directory = output_dir + "/" + genus
#     if not os.path.exists(directory):
#         os.makedirs(directory)
    
#     # Save the genome
#     filename = os.path.join(directory, f"{str(tax_id)}.fa")
#     with open(filename, "a") as output_file:
#         SeqIO.write(sequence, output_file, "fasta")


# Create a fasta file for the sequence
def save_to_file(sequence, output_dir, tax_id):

    # If the output directory doesn't exist, make it
    directory = output_dir
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Save the genome
    filename = os.path.join(directory, f"{str(tax_id)}.fa")
    with open(filename, "a") as output_file:
        SeqIO.write(sequence, output_file, "fasta")


# Main
if __name__ == "__main__":

    # Set parameters
    batch_size = 1000
    output_dir = "microbial_genomes/fasta_all_NC"
    input_file = "NCBI_sequenceIDs_NC.mapping"
    # genus_filter = ["Streptococcus",
    #                 "Mycobacterium",
    #                 "Staphylococcus"]
    
    # Require an email address argument
    parser = argparse.ArgumentParser(description="Download genomes from Refseq, filtered by genus.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    init_Entrez(args.email)
    
    # Keep track of how many of each genome we find
    # counts = {genus: 0 for genus in genus_filter}

    # Read in all ids
    with open(input_file, 'r') as file:
        refseq_ids = file.readlines()
    for i in range(len(refseq_ids)):
        refseq_ids[i] = refseq_ids[i].strip().split()

    # Filter out everything except NC_ ids
    # ID number 94200 is where I had the issue when using the data without the filter.
    # nc_ids = []
    # for id in refseq_ids:
    #     if "NC_" in id:
    #         nc_ids.append(id)
    # refseq_ids = nc_ids

    # Separate the list into batches
    num_ids = len(refseq_ids)
    for i in range(0, num_ids, batch_size):
        time.sleep(1)
        sys.stdout.write(f"\rProcessing sequence {i}/{num_ids}")
        batch = refseq_ids[i:i+batch_size]
        id_batch = [tuple[0] for tuple in batch]

        # Get the batch of sequences and summaries
        sequences = get_file_batch(id_batch)

        # summaries = get_summary_batch(id_batch)
        # # Create a list of tuples containing the accession number and the taxid
        # id_tuples = [(summary["AccessionVersion"], int(summary["TaxId"])) for summary in summaries]

        for id_tuple in batch:
            sequence = next(temp for temp in sequences if temp.id == id_tuple[0])
            # genus = sequence.description.split()[1]
            tax_id = id_tuple[1]

            save_to_file(sequence, output_dir, tax_id)

            # If it's a genus we want, save it and add 1 to its count
            # if genus in genus_filter:
                # save_to_file(sequence, output_dir, genus, tax_id)
                # counts[genus] += 1
    
    # Print results summary
    # sys.stdout.write(f"\rProcessed {num_ids} ids, saved:\n")
    # for genus in counts:
    #     print(f"{counts[genus]} {genus} sequences")
