# genomeWars
Repo for Bio 465 capstone project to test for genome contamination in microbial genomes. 

# Team Members 
    - Eliza Homer 
    - Kylee Bates 
    - Adam Rice 
    - Alex Brown 

# Instructions 

## Step 0: Downloading genomes 

1. Create a conda env with the required packages using the data_requirements.txt that is found in the data folder.
`conda create --name download_genomes --file data_requirements.txt`
`conda activate download_genomes`

2. Download the data for the 9 genera 
Sequence ids where found using this website: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=674962

Go into the `step_0_data` folder. Run the download_accession_number_9_genera.py script as follows: 

`python download_accession_numbers_9_genera.py --email your.email@email.com`

This script should take less than 10 minutes to download the genomes and it will download about 45 genomes. The accession numbers for this genomes can be found in the `accession_numbers_9_genera.csv` file. Alternatively, one can use the genomes that have already been downloaded in the `step_0_data/microbial_genomes_9_genera` folder. 

3. Download the data for the 17868 genomes

Go into the `step_0_data` folder. Run the download_microbe_genomes_fasta.py script to download the 17000 genomes to the folder called `microbial_genomes_17000`. This script takes about 5 hours to run. This script downloads the accession numbers found in the `NCBI_sequenceIDs_NC.txt` file. Note for TA: The genomes are already downloading into the `microbial_genomes_17000` folder on the supercomputer. 

`python download_microbe_genomes_fasta.py --email your.email@email.com`


## Step 1: Run the GUNC algorithm 

1. Set up a virtual environment with the required packages using the gunc_requirements.txt that is found in the GUNC folder. 
```sh
conda create --name gunc_env --file gunc_requirements.txt
conda activate gunc_env
```

- Clone the git repo for GUNC
    https://github.com/grp-bork/gunc

- Follow the documentation on the gunc website to set up the gunc algorithm 
    https://grp-bork.embl-community.io/gunc/

- One of those steps will be downloading the GUNC Database
    `gunc download_db GUNC_DB` 

- You will also need to install diamond 2.0.4. Follow the diamond installation documentation. 
    https://gensoft.pasteur.fr/docs/diamond/2.0.4/2_Installation.html

- You will also need to install prodigal. Do so through conda and/or follow the instructions on their website 
    `conda install bioconda::prodigal`
    https://anaconda.org/bioconda/prodigal

2. Run the GUNC algorithm for the 9 genera. 

Run the GUNC algorithm on the 9 genera by going into the `step_1_GUNC/run_output_9_genera` folder and submitting the multithreading.sh script. This script takes about 20 hours to run.

`sbatch submit_multithreading.sh`

After the slurm job had finsihed running, the `slurm_files` folder will contain the slurm files from the job, and the `output` folder will contain the raw output files. To extract the data from the raw files, run the `extract_data.py` script. This will go through all the raw output files and collect the data into the output.tsv file. The output.tsv will be used for further downstream anaylsis. There is an `expected_output_9_genera.tsv` file that contains the expected output from this step. 

```sh
python extract_data.py
```

3. Run the GUNC algorithm for the 17000 genomes. 

Run the GUNC algorithm on the 17000 genomes by going into the `step_1_GUNC/run_output_17000` folder and submitting the submit_multithreading.sh script. THis script takes about 150 hours to run.

`sbatch submit_multithreading.sh`

After the slurm job had finsihed running, the `slurm_files` folder will contain the slurm files from the job, and the `output` folder will contain the raw output files. To extract the data from the raw files, run the `extract_data.py` script. This will go through all the raw output files and collect the data into the output.tsv file. The output.tsv will be used for further downstream anaylsis. There is an `expected_output_17000_genera.tsv` file that contains the expected output from this step. 

```sh
python extract_data.py
```

## Step 2: Run VecScreen

### Step 2.1: Setup

1. Change your directory to `step_2_VecScreen`.

```sh
cd $HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen
```

2. Download BLAST+

Having the BLAST+ command line tool is necessary for running both the VecScreen and BLAST scripts. If you have cloned our repo, the executables for BLAST+ 2.15.0 should already be downloaded in a folder called `blast_executables` within this subfolder. The following command was run to obtain these files. If you wish to use a more recent version of BLAST, you will need to update the URL with the updated version which can be found on NCBI's website - https://www.ncbi.nlm.nih.gov/books/NBK52640/.

```sh
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz 
tar -xzf ncbi-blast-2.15.0+-x64-linux.tar.gz
mkdir blast_executables
mv ncbi-blast-2.15.0+* blast_executables/
```

3. Set path for BLAST+

```sh
export PATH=$PATH:$HOME/fsl_groups/grp_Bio465_GenomeWars/compute/genomeWars/step_2_VecScreen/blast_executables/ncbi-blast-2.15.0+/bin
```

4. Create UniVec DB

VecScreen uses blastn with certain parameters and aligns against a particular database, UniVec. First, we have to create the UniVec database to blast against. Again, if you have cloned this repo, you should not need to actually run these commands because the UniVec database should already be downloaded in a folder within this subfolder called `UniVecDB`. However, running the following commands will download the UniVec files and make them into a BLAST database.

```sh
wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
makeblastdb -in UniVec -dbtype nucl
mkdir UniVecDB
mv UniVec.* UniVecDB/
```

### Step 2.2: Run VecScreen

The following command is what is used in the scripts to run VecScreen with blastn.

```sh
blastn -db $db_path -query $file_path -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 5 -out $out_path
```

The following scripts will run VecScreen on the microbial genomes. The first script runs VecScreen on the 17,868 complete genomes. This script takes approximately 15 1/2 hours to run. The expected xml output for each genome is found in the `run_17000/expected_xml_results` folder. If you actually run the script with the following command, your results will be found in a folder called `run_17000/xml_results`.

```sh
sbatch run_17000/vecscreen_for_17000.sh
```

The second script runs VecScreen on 47 genomes that come from the 9 different genera that Salzberg noted as particuarly problematic. This script takes approximately 2 1/2 minutes to run. The expected xml output for each genome is found in the `run_9_genera/expected_xml_results` folder. If you actually run the script with the following command, your results will be found in a folder called `run_9_genera/xml_results`.

```sh
sbatch run_9_genera/vecscreen_for_9_genera.sh
```

### Step 2.3: Parse the VecScreen Output

The following scripts parse all of the xml files and place all relevant information into one csv file. To run these scripts, a conda environment must be created and activated with the BioPython package installed. The following commands will create an environment called `blast`, activate the environment, and install the package.

```sh
conda create --name blast_env
conda activate blast_env
conda install -c conda-forge biopython
```

Running the following python script will parse the `run_17000/expected_xml_results` xml files and out put a file called `run_17000/vecscreen_17000_output.tsv`. The expected output is found in a file called `run_17000/expected_vecscreen_17000_output.tsv`. Running this should take about 6 minutes.

```sh
python run_17000/parse_17000_xml.py
```

Running the following python script will parse the `run_9_genera/expected_xml_results` xml files and out put a file called `run_9_genera/vecscreen_9_genera_output.tsv`. The expected output is found in a file called `run_9_genera/expected_vecscreen_9_genera_output.tsv`. Running this should take less than 10 seconds.

```sh
python run_9_genera/parse_9_genera_xml.py
```

## Step 3: BLAST against HomoSapien

Go into the `step_3_BlastHomoSapien` folder for this step
```sh
cd step_3_BlastHomoSapien
```

### Step 3.1: Obtain the Human Genome
**\[TAs\]** Step 3.1 is the only part of step 3 where the data is not already completed for you. You'll have to wait a maximum of about 30 seconds for some of these commands to run. 

1. Obtain the zipped GRCh38.p14 Human Genome Assembly (GCF_000001405.40) from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. Multiple methods can be used to do this, but here we'll use curl. Run the following commands:
```sh
mkdir db
cd db
# ~20 seconds
curl https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz > GRCh38_latest_genomic.fna.gz
```

2. Unzip the file using gzip (it'll take a minute or two)
```sh
# ~30 seconds
gzip -d GRCh38_latest_genomic.fna.gz
```
You can make sure it worked by running `head -10 GRCh38_latest_genomic.fna` to show the first 10 lines. The first line will say ">NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly" followed by lines of "N's"

To see all the different sequences the FASTA file has, you can look for all the '>' symboles in the file by running `grep -F ">" GRCh38_latest_genomic.fna`

3. You should now have a file called "GRCh38_latest_genomic.fna" that is located in the db folder you created earlier. This file is about 3GB large. To make it into a database compatible with BLAST, use makeblastdb. 

```sh
# ~30 seconds
makeblastdb -in ./GRCh38_latest_genomic.fna -dbtype nucl -title GRCh38_latest_genomic -parse_seqids
```
*Note:* Before running this command, make sure you have access to the BLAST executables, and have exported it into your $PATH environment variable during the vecscreen set up. See step 2.1 number 3. 

This creates a number of files, including a .njs file containing the database metadata. You're now ready to run BLAST analysis against the human genome.

4. After this step, make sure to go back to the `step_3_BlastHomoSapien` folder. All commands remaining for step 3 will be done from this folder.
```sh
cd ..
```

### Step 3.2: Chunk the microbial data
Before finding alignments against the Human genome, we must chunk the query genomes (split them up into many different consecutive sequences of a specified length).

1. To chunk the genomes created in Step 0, run the chunk_folder.sh script using the following syntax:
* **\[TAs\]** This step takes a long time, and has already been completed. The output folders are 

`step_3_BlastHomoSapien/pre_chunked_9_genera` and `step_3_BlastHomoSapien/pre_chunked_17000`
```sh
# ~10 seconds
sh chunk_folder.sh ../step_0_data/microbial_genomes_9_genera ./chunked_9_genera
# ~2 hours
sh chunk_folder.sh ../step_0_data/microbial_genomes_17000 ./chunked_17000
```

### Step 3.3: Run BLAST Alignment Algorithm
To run the BLAST alignement algortihm against the human genome database we created, run the following command:
* **\[TAs\]** This step takes a long time, and has already been completed. The output folders are `step_3_BlastHomoSapien/pre_results_9_genera` and `step_3_BlastHomoSapien/pre_results_17000`
```sh
# ~10 minutes
sh human_blastn.sh -q ./pre_chunked_9_genera -r ./results_9_genera
# ~ 25 hours
sh human_blastn.sh -q ./pre_chunked_17000 -r ./results_17000
```

### Step 3.4: Parse the BLAST Output
Once the data has been processed into XML files, we'll parse it using a python script like this:
* This script requires that you use the Bio module, which you should have imported into a conda blast-env during step 2
```sh
# ~5 seconds
python3 parse_xml_output.py ./pre_results_9_genera ./_parsed/parsed_9_genera.tsv
# ~15 minutes
python3 parse_xml_output.py ./pre_results_17000 ./_parsed/parsed_17000.tsv
```
* **\[TAs\]** This step takes a long time, and has already been completed. The output paths are `step_3_BlastHomoSapien/_parsed/pre_parsed_9_genera.tsv` and `step_3_BlastHomoSapien/_parsed/pre_parsed_17000.tsv`


## Step 4: Run FCS-GX


## Step 5: Create Figures

1. Go to the `step_5_figures` directory. 

```sh
cd step_5_figures
```

2. Create figure 1

Go to the `step_5_figures/figure_1_GUNC` folder.

```sh
cd step_5_figures/figure_1_GUNC
```

The `generate_figure_1.R` script creates Figure 1, which corresponds to the GUNC output for the 9 genera. Download the `step_1_GUNC/run_output_17000/expected_output_17000.tsv` file and run the script in RStudio. Make sure your working directory in RStudio is set to the folder where you have placed the output file. The expected output for this figure is found in the `step_5_figures/figure_1_GUNC/` folder as `expected_figure_1.png`.

This figures should the names of the genomes that had 25% or more contamination using the GUNC algorithm. 

3. Create figure 2

Go to the `step_5_figures/figure_2_vecscreen` folder.

```sh
cd step_5_figures/figure_2_vecscreen
```

The `generate_figure_2.R` script creates Figure 2, which corresponds to the VecScreen output. Download the `step_2_VecScreen/run_17000/expected_vecscreen_17000_output.tsv` file and run the script in RStudio. Make sure your working directory in RStudio is set to the folder where you have placed the output file. The expected output for this figure is found in the `step_5_figures/figure_2_vecscreen` folder as `expected_figure_2.png`.

4. Create figure 3 

Go to the `step_5_figures/figure_3_blast` folder.

```sh
cd step_5_figures/figure_3_blast
```

The `generate_figure_3.R` script creates  Figure 3, which corresponds to the BLAST output. Download the `step_3_BlastHomoSapien/_parsed/pre_parsed_17000.tsv` file and run the script in RStudio. Make sure your working directory in RStudio is set to the folder where you have placed the output file. The expected output for this figure is found in the `step_5_figures/figure_3_blast` folder as `expected_figure_3.png`.


5. Create figures 4a and 4b

Go to the `step_5_figures/figure_4_FCS-GX` folder.

```sh
cd step_5_figures/figure_4_FCS-GX
```

The `generate_figure_4.R` script creates figures 4a and 4b, which corresponds to the FCS-GX output. Download the `step_5_figures/figure_4_FCS-GX/summary.txt` file and run the script in RStudio. Make sure your working directory in RStudio is set to the folder where you have placed the output file. The expected output for these figures are found in the `step_5_figures/figure_4_FCS-GX` folder as `expected_figure_4a.png` and `expected_figure_4a.png`.