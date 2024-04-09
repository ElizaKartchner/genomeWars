#!/bin/bash

#SBATCH --time=20:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --cpus-per-task=8     # number of CPU cores per task
#SBATCH --nodes=4   # number of nodes
#SBATCH --mem=300G             # total memory for the job
#SBATCH -J "9_gunc_run"   # job name
#SBATCH --mail-user=ek252@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output "slurm_files/output_.%A.out"   # name of stdout output file
#***** sg grp_Bio465_GenomeWars 'sbatch submit_multithreading.sh'

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
echo "Loading parallel"
module load parallel

export NUM_TASKS=16
export NUM_THREADS=8

# NOTE: verify that this is the correct input dir
export GENOMES_DIRECTORY="../../step_0_data/microbial_genomes_9_genera"

export TSV_FILE="GUNC.progenomes_2.1.maxCSS_level.tsv"
export NUMBER_FILE_NOT_EXIST=0

# NOTE: verify that this is the correct destination folder
export DESTINATION_FOLDER="output"
export DESTINATION_FILE="output.tsv"

echo "Starting..."

run_gunc() {
    echo "starting gunc func"
    local file="$1"

    filename=$(basename "$file")
    fa_name="${filename%.*}"
    out_dir="${DESTINATION_FOLDER}/output_${fa_name}"

    [ ! -d $out_dir ] && mkdir -p $out_dir

    echo "Screening $fa_name..."

    gunc run -i "$file" -r ../gunc_db_progenomes2.1.dmnd --threads 8 --out_dir $out_dir

}
export -f run_gunc

echo "beginning parallel..."
parallel --jobs $NUM_THREADS --keep-order run_gunc {} ::: "$GENOMES_DIRECTORY"/*

echo "Finished"

