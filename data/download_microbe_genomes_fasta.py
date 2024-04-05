""" This file is used to download about 18,000 genomes. It takes about 4.5 hours to do so."""

import os
import argparse
from Bio import Entrez, SeqIO 


from Bio import Entrez, SeqIO

def download_genomes(email, refseq_id):
    # Set your email address for Entrez
    Entrez.email = email

    # Use Entrez to search for the specified term and retrieve the list of IDs
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="fasta", retmode="text")
    sequence = SeqIO.read(handle, "fasta")
    handle.close()

    # Download and save each genome
    output_dir = f"microbial_genomes_17000"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filename = os.path.join(output_dir, f"{sequence.id.replace(".", "_")}.fa")
    with open(filename, "w") as output_file:
        SeqIO.write(sequence, output_file, "fasta")
    print(f"Genome {sequence.id} downloaded and saved as {filename}")


if __name__ == "__main__":
    # Set your email address 
    parser = argparse.ArgumentParser(description="Download genomes from Refseq.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    email_address = args.email
    
    refseq_ids = []
    input_file = "NCBI_sequenceIDs_NC.txt"
    with open(input_file, "r") as file:
        refseq_ids = file.read().splitlines()

    print(refseq_ids[1:10])
    print(len(refseq_ids))

    for refseq_id in refseq_ids:
        download_genomes(email_address, refseq_id)
