#!/bin/bash

clear_line() {
  printf "\r\033[K"
}



data_folder=$1 # "../../Step1_Acquire_Data/microbial_genomes"
chunked_data_folder=$2 # "../Step2.4/query"
chunk_size=1000

mkdir $chunked_data_folder 1>>/dev/null

count=0
total=$(ls -l "$data_folder" | grep "^-" | wc -l)

echo "Chunking files"
printf "..."

for file_path in "$data_folder"/*; do
  if [ -f "$file_path" ]; then
    clear_line
    printf "$count/$total"
    file_name=$(basename "$file_path")
    out_path="$chunked_data_folder/$file_name.chk.fa"

    python3 chunk_genomes.py $file_path $out_path $chunk_size 1>>/dev/null

    # blastn -db $db_path -query $file_path -out $out_path
    # grep -F ">" "$out_path"
    ((count++))
  fi
done

clear_line
echo "chunked $count files"