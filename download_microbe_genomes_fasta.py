import os
import argparse
from Bio import Entrez, SeqIO 

# Leucothrix, Thalassomonas, Velarivirus, Tritimovirus, Dinovernavirus, Bacillarnavirus, Rymovirus, Ignicoccus, Salinimicrobium

# NC_014122.1

from Bio import Entrez, SeqIO

def download_genomes(email, refseq_id):
    # Set your email address for Entrez
    Entrez.email = email

    # Use Entrez to search for the specified term and retrieve the list of IDs
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="fasta", retmode="text")
    sequence = SeqIO.read(handle, "fasta")
    handle.close()


    # Download and save each genome
    output_dir = f"microbial_genomes/fasta_files"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filename = os.path.join(output_dir, f"{sequence.id}.fasta")
    with open(filename, "w") as output_file:
        SeqIO.write(sequence, output_file, "fasta")
    print(f"Genome {sequence.id} downloaded and saved as {filename}")


if __name__ == "__main__":
    # Set your email address and search term
    parser = argparse.ArgumentParser(description="Download genomes from Refseq.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    email_address = args.email
    
    refseq_ids = ["NC_014122.1"]

    for refeq_id in refseq_ids:
        download_genomes(email_address, refeq_id)
