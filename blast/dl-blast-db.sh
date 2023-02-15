#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-blast-db-%j.out

# Command-line args
target_dir=$1

update_blastdb_path=/home/jelmer/miniconda3/envs/blast-env/bin/update_blastdb.pl

# Bash strict mode
set -euo pipefail

# Report
echo "## Starting script blast-db.sh..."
echo "## This script will download a Blast database."
date
echo "## Command-line args:"
echo "## Input FASTA file: $target_dir"
echo "## Path to NCBI's update_blast_db.pl script: $update_blastdb_path"
echo -e "---------------\n\n"

# Download Blast DB:
mkdir -p "$target_dir" 
cd "$target_dir" || exit

"$update_blastdb_path" --decompress --num_threads "$SLURM_CPUS_ON_NODE" nt
"$update_blastdb_path" --decompress taxdb

# Report
echo "## Done with script blast-db.sh"
date