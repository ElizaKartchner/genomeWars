#!/bin/bash

#SBATCH --time=0:10:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "run_9_genera"   # job name

_p=" >"

db_path="$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/UniVecDB/UniVec"
query_folder="$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_0_data/microbial_genomes_9_genera/"
results_folder="$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/run_9_genera/xml_results"

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

echo "DB: $db_path"
echo "Query: $query_folder"
echo "Results: $results_folder"

# export BLAST=$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/BLAST/blast_executables/ncbi-blast-2.15.0+/bin

# Check if the folder already exists
if [ -d "$results_folder" ]; then
  rm -r "$results_folder"
fi

# Create the folder
mkdir "$results_folder"
echo "$_p Results Folder: $results_folder."

highlight='\033[0;31m'
normal='\033[0m'
count=0

file_count=$(ls -l "$query_folder" | grep "^-" | wc -l)

for file_path in "$query_folder"/*; do
  if [ -f "$file_path" ]; then
    file_name=$(basename "$file_path")
    out_path="$results_folder/$file_name.xml"
    ((count++))
    echo -e "${_p}${highlight} Processing file ($count of $file_count): $file_name $normal"

    blastn -db $db_path -query $file_path -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 5 -out $out_path
    grep -F ">" "$out_path"
  fi
done

echo "processed $count files"