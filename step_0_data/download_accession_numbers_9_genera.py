import os
import argparse
import pandas as pd
from Bio import Entrez, SeqIO 

from Bio import Entrez, SeqIO

def download_refseq_genomes(email, refseq_id):
    # Set your email address for Entrez
    Entrez.email = email

    # Use Entrez to search for the specified term and retrieve the list of IDs
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="fasta", retmode="text")
    sequence = SeqIO.read(handle, "fasta")
    handle.close()

    # Download and save each genome
    output_dir = f"microbial_genomes_9_genera"
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
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv("accession_numbers_9_genera.csv")

    unique_refseq_nums = df["refseq_num"].unique()
    unique_refseq_nums = unique_refseq_nums[unique_refseq_nums != 'na']

    print(unique_refseq_nums)
    print(len(unique_refseq_nums))

    for refseq_id in unique_refseq_nums:
        print(f"Now downloading refseq id: {refseq_id}")
        download_refseq_genomes(email_address, refseq_id)

