
clear_line() {
  printf "\r\033[K"
}

_p=" >"

db_path="db/GRCh38_latest_genomic.fna"
query_folder="query"
results_folder="results"

# Parse command line arguments
while getopts ":q:r:" opt; do
  case ${opt} in
    q )
      query_folder=$OPTARG
      ;;
    r )
      results_folder=$OPTARG
      ;;
    \? )
      echo "Usage: $0 [-q query] [-r results]"
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

summary_file="$results_folder/__summary.out"

printf "\n\n"
echo "DB: $db_path"
echo "Query: $query_folder"
echo "Results: $results_folder"
echo "Summary: $summary_file"
printf "\n"

# export PATH=$PATH:/home/awr28/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/BLAST/blast_executables/ncbi-blast-2.15.0+/bin

# Check if the folder already exists
if [ -d "$results_folder" ]; then
  rm -r "$results_folder"
fi

# Create the folder
mkdir "$results_folder"
count=0
alignment_count=0

echo "Alignments have been found in the following files" > $summary_file

file_count=$(ls -l "$query_folder" | grep "^-" | wc -l)

for file_path in "$query_folder"/*; do
  if [ -f "$file_path" ]; then
    file_name=$(basename "$file_path")
    out_path="$results_folder/$file_name.out.xml"
    clear_line
    printf "[$count/$file_count] [found: $alignment_count] [$file_name]"

    blastn -db $db_path -query $file_path -out $out_path -outfmt 5
    if grep -F "|" "$out_path" >>/dev/null; then
      ((alignment_count++))
      echo $file_name >> $summary_file
    fi
    ((count++))
  fi
done

clear_line
printf "\n\n"
echo "Processed $count files, found $alignment_count alignments"