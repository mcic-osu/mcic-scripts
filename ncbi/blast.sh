#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-blast-%j.out

## Help function
Help() {
    echo
    echo "$0: Run BLAST locally or remotely for a FASTA file."
    echo
    echo "Syntax: $0 -i <input-fasta> -o <output-file> ..."
    echo
    echo "Required options:"
    echo "   -i FILE      Input file ('query', a FASTA file)"
    echo "   -o FILE      Output file"
    echo
    echo "Other options:"
    echo "   -d STRING    Blast DB                              [default: 'nt']"
    echo "                If remote, e.g. 'nt' or 'nr'."
    echo "                If local, make sure to specify the the dir AND the DB name"
    echo "                e.g. '/fs/project/PAS0471/blast/nr-db/nr'"
    echo "   -h           Print this help message and exit"
    echo "   -l           Run BLAST locally                     [default: run BLAST remotely]"
    echo
    echo "Example:        $0 -i data/my.fa -o results/blast/blast.out -d /fs/project/PAS0471/blast/nt-db/nt -l"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}


# PARSE OPTIONS ----------------------------------------------------------------
## Option defaults
query_fa=""
blast_out=""
blast_db="nt"
remote="true"

## Parse command-line options
while getopts ':i:o:d:lh' flag; do
    case "${flag}" in
    i) query_fa="$OPTARG" ;;
    o) blast_out="$OPTARG" ;;
    d) blast_db="$OPTARG" ;;
    l) remote=false ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/blast-env
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08

## Bash strict mode
set -euo pipefail

## Process options
outdir=$(dirname "$blast_out")
mkdir -p "$outdir"

## Number of cores
n_cores=$SLURM_CPUS_PER_TASK

## Make paths absolute
[[ ! "$query_fa" =~ ^/ ]] && query_fa="$PWD"/"$query_fa"
[[ ! "$blast_out" =~ ^/ ]] && blast_out="$PWD"/"$blast_out"
[[ "$remote" = false && ! "$blast_db" =~ ^/ ]] && blast_db="$PWD"/"$blast_db" 

## Test input
[[ $query_fa = "" ]] && echo "## $0: ERROR: No input file (-i) provided" >&2 && exit 1
[[ $blast_out = "" ]] && echo "## $0: ERROR: No output file (-o) provided" >&2 && exit 1
[[ ! -f $query_fa ]] && echo "## $0: ERROR: Input file $query_fa does not exist" >&2 && exit 1

## Report
echo
echo "## Starting script blast.sh..."
date
echo
echo "## Input FASTA file:            $query_fa"
echo "## BLAST output file:           $blast_out"
echo "## BLAST database:              $blast_db"
echo "## Run blast remotely?          $remote"
[[ "$remote" = false ]] && echo "## Number of cores:             $n_cores"
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
        -num_threads "$n_cores"

fi


# WRAP UP ----------------------------------------------------------------------
echo -e "------------------------\n"
echo "## Listing the output file:"
ls -lh "$blast_out"
echo
echo "## Done with script blast.sh"
date
echo
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
