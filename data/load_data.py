from Bio import Entrez, SeqIO
import os

def download_refseq_by_genus(genus, output_folder):
    Entrez.email = "eliza.kartchner@gmail.com"  # Set your email address

    # Step 1: Search for the genus in the taxonomy database
    handle = Entrez.esearch(db="taxonomy", term=f"{genus}[Organism]", retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    if len(record["IdList"]) == 0:
        print(f"No records found for genus: {genus}")
        return
    else:
        print(f"{len(record["IdList"])} records found for genus: {genus}")

    
    refseq_records = []
    for tax_id in record["IdList"]:

        # Step 2: Fetch RefSeq nucleotide records for the given taxonomy ID
        #handle = Entrez.esearch(db="nucleotide", term=f"{tax_id}[TaxID] AND refseq[filter]", retmode="xml")  # This looks just as refseq
        handle = Entrez.esearch(db="nucleotide", term=f"{tax_id}[TaxID]", retmode="xml")  # This give us a LOT more files to look at 91
        record = Entrez.read(handle)
        handle.close()
        refseq_records.append(record)


    if len(refseq_records) == 0:
        print(f"No RefSeq records found for genus: {genus}")
        return
    else:
        print(f"{len(refseq_records)} RefSeq records found for genus: {genus}")

    # Step 3: Download and save each sequence in a separate file
    for record in refseq_records:
        for seq_id in record["IdList"]:
            with Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text") as handle:
                seq_record = SeqIO.read(handle, "fasta")

            output_filename = os.path.join(output_folder, f"{seq_record.id}.fa")

            with open(output_filename, "w") as output_file:
                SeqIO.write(seq_record, output_file, "fasta")

            print(f"Saved {seq_record.id} to {output_filename}")

if __name__ == "__main__":
    genus_name = "bacillarnavirus"   # Replace with your desired genus name
    output_directory = f"{genus_name}"  

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    download_refseq_by_genus(genus_name, output_directory)

