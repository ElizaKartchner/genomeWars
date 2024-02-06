import os
import argparse
from Bio import Entrez, SeqIO 

# Leucothrix, Thalassomonas, Velarivirus, Tritimovirus, Dinovernavirus, Bacillarnavirus, Rymovirus, Ignicoccus, Salinimicrobium

from Bio import Entrez, SeqIO

def download_genomes(email, microbe, num_results):
    # Set your email address for Entrez
    Entrez.email = email

    search_term = microbe + "[Organism]"

    # Use Entrez to search for the specified term and retrieve the list of IDs
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=num_results, idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    # Download and save each genome
    for uid in record["IdList"]:
        # Fetch the record
        handle = Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()

        # Save the record to a file
        output_dir = f"microbial_genomes/{microbe}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        filename = os.path.join(output_dir, f"{seq_record.id}.gbk")
        with open(filename, "w") as output_file:
            SeqIO.write(seq_record, output_file, "genbank")
        print(f"Genome {seq_record.id} downloaded and saved as {filename}")

if __name__ == "__main__":
    # Set your email address and search term
    parser = argparse.ArgumentParser(description="Download genomes from Refseq.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    email_address = args.email
    
    num_results = 1  # Number of genomes to download

    top_microbes = ["leucothrix", "thalassomonas", "velarivirus", "tritimovirus", "dinovernavirus", "bacillarnavirus", "rymovirus", "ignicoccus", "salinimicrobium"]

    for microbe in top_microbes:
        download_genomes(email_address, microbe, num_results)
