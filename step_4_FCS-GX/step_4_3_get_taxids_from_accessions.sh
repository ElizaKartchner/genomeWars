#!/bin/bash

# Check usage
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_ids_file> <output_map>"
    echo "  - input_ids_file: Text file containing accession IDs"
    echo "  - output_map: File which will contain the map between filename, accession, and taxid"
    exit 1
fi

# Get list of accession numbers from the command line
export SEQ_IDS="$1"
export OUT_MAP="$2"

# Define the batch size
batch_size=100

# Create a temporary subfolder
export TEMP_DIR="temp"
mkdir -p $TEMP_DIR

# Split the input file into chunks in the temporary subfolder
export CHUNK_DIR=$TEMP_DIR/seq_id_chunks
mkdir -p $CHUNK_DIR
split -l $batch_size $SEQ_IDS $CHUNK_DIR/

# Get the total number of chunks
total_chunks=$(ls $CHUNK_DIR/* | wc -l)

# Process each chunk
count=0
for file in $CHUNK_DIR/*
do
    # Take the file with filenames and accessions and create a file with accessions and taxids
    awk '{print $2}' "$file" > temp1.temp
    cat temp1.temp | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element AccessionVersion,TaxId > temp2.temp
    
    # Perform an inner join on the taxid column to create one large mapping file
    sort -k2 $file > temp3.temp
    sort -k1 temp2.temp > temp4.temp
    join -1 2 -2 1 temp3.temp temp4.temp | awk '{print $2, $1, $3}' > "$file.temp"

    # Update the count
    count=$((count + 1))
    echo -ne "\rProcessed chunk $count of $total_chunks"
done
echo -ne "\r"
echo -ne "\033[K"
echo -e "Processed all $total_chunks chunks"
rm temp1.temp
rm temp2.temp
rm temp3.temp
rm temp4.temp

# Combine all chunk results into one file
cat $CHUNK_DIR/*.temp >> $OUT_MAP

# Delete the temporary subfolder and the chunks
rm -rf $CHUNK_DIR
rm -rf $TEMP_DIR

echo "Mapping file created"