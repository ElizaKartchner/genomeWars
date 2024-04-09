#!/bin/bash

# Set Project Directories
export PROJECT_DIR="$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_4_FCS-GX"
export VIRTUAL_DIR="/dev/shm"
export REPORTS_DIR="$PROJECT_DIR/reports"
export LOGGING_DIR="$PROJECT_DIR/logs/gx"
export SLURM_DIR="$PROJECT_DIR/logs/slurm"
export SUMMARY_CONDENSED_FILE="$PROJECT_DIR/summary.txt"
export SUMMARY_FULL_FILE="$PROJECT_DIR/summary.rpt"

# Default supercomputer parameters
export CORES_PER_NODE=32
export THREADS_PER_NODE=8
export CORES_PER_THREAD=4

# Parse command-line options
while getopts ":c:t:p:" opt; do
  case "$opt" in
    c) CORES_PER_NODE="$OPTARG" ;;
    t) THREADS_PER_NODE="$OPTARG" ;;
    p) CORES_PER_THREAD="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1 ;;
  esac
done
shift $((OPTIND-1))  # Shift past the options

# Check usage
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [-c cores] [-t threads] [-p cores_per_thread] [-o out_sum] <gen_dir> <taxid_map> <email>"
    echo "  -c cores: Number of cores per node (default: 1)"
    echo "  -t threads: Number of threads per node (default: 1)"
    echo "  -p cores_per_thread: Number of cores per thread (default: 1)"
    echo "  - gen_dir: The directory that contains the fasta files to screen"
    echo "  - taxid_map: The mapping file of filenames, accessions, and taxids"
    echo "  - email: The email of the user to run the slurm script"
    exit 1
fi

# Input parameters
export INPUT_DIR="$1"
export TAX_ID_FILE="$2"
export EMAIL_ADDRESS="$3"

# Create Project Directories
[ ! -d "$REPORTS_DIR" ] && mkdir -p "$REPORTS_DIR"
[ ! -d "$LOGGING_DIR" ] && mkdir -p "$LOGGING_DIR"
[ ! -d "$SLURM_DIR" ] && mkdir -p "$SLURM_DIR"

sbatch \
    --tasks $CORES_PER_NODE \
    --mail-user $EMAIL_ADDRESS \
    --output "$SLURM_DIR/FCS-GX.%A_%a.out" \
    "$PROJECT_DIR/core/screen_genomes.sh"
