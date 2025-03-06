#!/bin/bash

# Building a SNPeff database
# See https://pcingola.github.io/SnpEff/snpeff/build_db

# Load the software
conda_env_type=/fs/ess/PAS0471/jelmer/conda/snpeff
module load miniconda3
conda activate "$conda_env"

# Genome characteristics
genome_id=GCF_002878395.1
species_name="Capsicum annuum"

# Derived variables
data_dir="$conda_env"/share/snpeff-5.2-1/data
config="$conda_env"/share/snpeff-5.2-1/snpEff.config
species_name2=${species_name/ /_}                       # E.g. "Capsicum_annuum"

# Test
[[ ! -d "$data_dir" ]] && echo "ERROR: SnpEff data dir does not exist"
[[ ! -f "$config" ]] && echo "ERROR: SnpEff config file does not exist"

# Download the reference genome files
sbatch mcic-scripts/download/dl-genomes.sh -A "$genome_id" -o data/ref

# Copy the genome assembly and annotation into the SnpEff data dir
mkdir -p "$data_dir"/"$genome_id"
cp -v data/ref/"$genome_id".fna "$data_dir"/"$genome_id"/sequences.fa
cp -v data/ref/"$genome_id".gtf "$data_dir"/"$genome_id"/genes.gtf
cp -v data/ref/"$genome_id".faa "$data_dir"/"$genome_id"/protein.fa

# Edit the snpEff.config file
echo -e "\n# $species_name genome, version $genome_id" >> "$config"
echo "${genome_id}.genome : $species_name2" >> "$config"

# Build the SnpEff database
sbatch -A PAS0471 -c 4 --wrap="snpEff build -gtf22 -verbose $genome_id"
