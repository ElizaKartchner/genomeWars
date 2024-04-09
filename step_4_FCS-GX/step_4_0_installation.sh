#!/bin/bash

export PROJECT_DIR="$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_4_FCS-GX"

# Install Entrez Direct to get the taxids later
echo "Installing Entrez Direct"
yes | sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

# Download FCS-GX python wrapper and singularity file
echo "Installing FCS-GX"
mkdir -p "$PROJECT_DIR/core/full_database"
curl -L https://github.com/ncbi/fcs/raw/main/dist/fcs.py -o "$PROJECT_DIR/core/fcs.py"
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif -Lo "$PROJECT_DIR/core/fcs-gx.sif"
export FCS_DEFAULT_IMAGE=fcs-gx.sif
echo "Loading singularity"
module load singularity

# Parse the command-line options
if [ "$1" == "-f" ]; then
    echo "Downloading full FCS-GX database..."
    cd "$PROJECT_DIR/core"
    SOURCE_DB_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest"
    LOCAL_DB="$PROJECT_DIR/core/full_database"
    python3 $PROJECT_DIR/core/fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb" 
    cd -
    echo "Database downloaded"
fi