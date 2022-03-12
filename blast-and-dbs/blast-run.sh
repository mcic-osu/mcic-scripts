#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-blast-%j.out

# SETUP ------------------------------------------------------------------------
## Load software
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate blast-env

## Bash strict mode
set -euo pipefail

## Help
Help() {
    # Display Help
    echo
    echo "## $0: Run BLAST for an input (query) FASTA file, either locally or remotely."
    echo
    echo "## Syntax: $0 -i <input-fasta> -o <output> -d <blast-db> [-r] [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     Input file ('query', a FASTA file) (REQUIRED)"
    echo "## -o     Output file (REQUIRED)"
    echo "## -d     Blast DB. If remote, e.g. 'nt' or 'nr'. If local, make sure to specify the dir AND THE DB name, e.g. refdata/blast/nr-db/nr"
    echo "## -r     Run BLAST remotely (default: run BLAST locally)"
    echo "## Example: $0 -q my.fasta -o results/blast/blast.out -d refdata/blast/nt-db/nt"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
query_fa=""      # BLAST input file name (FASTA)
blast_out=""     # BLAST output file name
blast_db="nr"
remote="false"

## Parse command-line options
while getopts ':i:o:d:rh' flag; do
    case "${flag}" in
    i) query_fa="$OPTARG" ;;
    o) blast_out="$OPTARG" ;;
    d) blast_db="$OPTARG" ;;
    r) remote=true ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process options
outdir=$(dirname "$blast_out")
mkdir -p "$outdir"

## Test input
[[ $query_fa = "" ]] && echo "## $0: ERROR: No input file (-i) provided" >&2 && exit 1
[[ $blast_out = "" ]] && echo "## $0: ERROR: No output file (-o) provided" >&2 && exit 1
[[ ! -f $query_fa ]] && echo "## $0: ERROR: Input file $query_fa does not exist" >&2 && exit 1

## Report
echo "## Starting script blast-run..."
date
echo "## Command-line args:"
echo "## Input FASTA file:     $query_fa"
echo "## BLAST output file:    $blast_out"
echo "## BLAST database:       $blast_db"
echo "## Run blast remotely?   $remote"
echo -e "--------------------\n"


# RUN BLAST --------------------------------------------------------------------
if [ "$remote" != false ]; then

    echo "## Running BLAST *remotely*..."
    blastn \
        -remote \
        -task blastn \
        -db "$blast_db" \
        -query "$query_fa" \
        -out "$blast_out" \
        -outfmt 6 -evalue 1e-6

else

    echo "## Running BLAST *locally*..."

    blast_db_dir=$(dirname "$blast_db")
    blast_db_name=$(basename "$blast_db")

    [[ ! -d "$blast_db_dir" ]] && echo "## $0: ERROR: Blast DB dir $blast_db_dir does not exist" >&2 && exit 1

    cd "$blast_db_dir"

    blastn \
        -task blastn \
        -db "$blast_db_name" \
        -query "$query_fa" \
        -out "$blast_out" \
        -outfmt 6 -evalue 1e-6 \
        -num_threads "$SLURM_CPUS_ON_NODE"

fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n\n## Listing the output file:"
ls -lh "$blast_out"

echo -e "\n## Done with script blast-run.sh"
date
