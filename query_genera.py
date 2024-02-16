import os
import argparse
from Bio import Entrez, SeqIO 

# Leucothrix, Thalassomonas, Velarivirus, Tritimovirus, Dinovernavirus, Bacillarnavirus, Rymovirus, Ignicoccus, Salinimicrobium

# NC_014122.1

from Bio import Entrez, SeqIO


def init_Entrez(email):
    # Set your email address for Entrez
    Entrez.email = email

def get_batch(refseq_ids):
    # Use Entrez to search for the specified term and retrieve the list of IDs
    handle = Entrez.efetch(db="nucleotide", id=refseq_ids, rettype="fasta", retmode="text")
    sequences = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return sequences

def save_to_file(sequence, output_dir):
    # If the output directory doesn't exist, make it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Download and save the genome
    filename = os.path.join(output_dir, f"{sequence.id.replace(".", "_")}.fa")
    with open(filename, "w") as output_file:
        SeqIO.write(sequence, output_file, "fasta")
    print(f"Genome {sequence.id} downloaded and saved as {filename}")


if __name__ == "__main__":
    # Require an email address argument
    parser = argparse.ArgumentParser(description="Download genomes from Refseq, filtered by genus.")
    parser.add_argument("--email", required=True, help="Your email address for Entrez.")
    args = parser.parse_args()
    init_Entrez(args.email)

    output_dir = "microbial_genomes/filter_test"
    input_file = "NCBI_sequenceIDs.txt"
    genus_filter = ["Streptococcus",
                    "Mycobacterium",
                    "Staphylococcus"]
    
    with open(input_file, 'r') as file:
        refseq_ids = file.readlines()

    # refseq_ids = ["NC_014122.1",
    #               "NZ_CP019964.1",
    #               "NZ_CP030847.1",
    #               "NZ_AP011528.1",
    #               "NZ_CP023154.1",
    #               "NC_013741.1",
    #               "NC_008213.1",
    #               "NC_021058.1",
    #               "NZ_CP008888.1",
    #               "NC_012632.1"]

    filtered_refseq_ids = []

    
    sequences = get_batch(refseq_ids)

    for sequence in sequences:
        genus = sequence.description.split()[1]

        if genus in genus_filter:
            filtered_refseq_ids.append(sequence.id)
            save_to_file(sequence, output_dir)

    print("Filtered IDs:")
    for id in filtered_refseq_ids:
        print(id)
     