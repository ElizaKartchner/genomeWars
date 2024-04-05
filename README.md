# genomeWars
Repo for Bio 465 capstone project to test for genome contamination in microbial genomes. 

# Team Member 
    - Eliza Homer 
    - Kylee Bates 
    - Adam Rice 
    - Alex Brown 

# Instructions 

## Downloading genomes 

1. Create a conda env with the required packages using the data_requirements.txt that is found in the data folder.
`conda create --name download_genomes --file data_requirements.txt`
`conda activate download_genomes`

2. Download the data for the 9 genera 
Sequence ids where found using this website: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=674962

Go into the data folder. Run the download_accession_number_9_genera.py script as follows: 

`python download_accession_numbers_9_genera.py --email your.email@email.com`

This script should take less than 10 minutes to download the genomes and it will download about 45 genomes. Alternatively, one can use the genomes that have already been downloaded in the `data/microbial_genomes_9_genera` folder. 

3. Download the data for the 17868 genomes

Go into the data folder. Run the download_microbe_genomes_fasta.py script to download the 17000 genomes. 

`python download_microbe_genomes_fasta.py --email your.email@email.com`


## Run the GUNC algorithm 

1. Set up a virtual environment with the required packages using the gunc_requirements.txt that is found in the GUNC folder. 
`conda create --name gunc_env --file gunc_requirements.txt`
`conda activate gunc_env`

- Clone the git repo for GUNC
    https://github.com/grp-bork/gunc

- Follow the documentation on the gunc website to set up the gunc algorithm 
    https://grp-bork.embl-community.io/gunc/

- One of those steps will be downloading the GUNC Database
    gunc download_db GUNC_DB 

- You will also need to install diamond 2.0.4. Follow the diamond installation documentation. 
    https://gensoft.pasteur.fr/docs/diamond/2.0.4/2_Installation.html

- You will also need to install prodigal. Do so through conda and/or follow the instructions on thier website 
    conda install bioconda::prodigal
    https://anaconda.org/bioconda/prodigal




