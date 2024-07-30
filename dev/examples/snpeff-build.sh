#!/bin/bash

# ==============================================================================
#                    2024-06-23 -- Building a SNPeff database
# ==============================================================================
# /fs/ess/PAS0471/jelmer/assist/00_misc/2023-12_chee

#? See https://pcingola.github.io/SnpEff/snpeff/build_db
conda_env=/fs/ess/PAS0471/jelmer/conda/snpeff
genome_id=GCF_002878395.1

# Download the reference genome files (you can do this any which way, but with this script is the easiest for me)
sbatch mcic-scripts/download/dl-genomes.sh -A "$genome_id" -o data/ref

# Copy the genome assembly and annotation into the SnpEff data dir
data_dir="$conda_env"/share/snpeff-5.2-1/data
mkdir -p "$data_dir"/"$genome_id"
cp -v data/ref/"$genome_id".fna "$data_dir"/"$genome_id"/sequences.fa
cp -v data/ref/"$genome_id".gtf "$data_dir"/"$genome_id"/genes.gtf
cp -v data/ref/"$genome_id".faa "$data_dir"/"$genome_id"/protein.fa

# Edit the snpEff.config file
config="$conda_env"/share/snpeff-5.2-1/snpEff.config
echo -e "\n# Capsicum annuum genome, version $genome_id" >> "$config"
echo "$genome_id : Capsicum_annuum" >> "$config"

# Build the SnpEff database
conda activate "$conda_env"
sbatch -A PAS0471 -c 4 --wrap="snpEff build -gtf22 -verbose GCF_002878395.1"
