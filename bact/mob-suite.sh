#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=mobsuite
#SBATCH --output=slurm-mobsuite-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run mob-suite to detect plasmids in a genome."
    echo
    echo "Syntax: $0 -i <genome-FASTA> -o <output-dir> ..."
    echo
    echo "Required options:"
    echo "-i FILE           Genome (nucleotide) FASTA file"
    echo "-o DIR            Output directory"
    echo
    echo "Other options:"
    echo "-h                Print this help message and exit"
    echo
    echo "Example:          $0 -i my_genome.fa -o results/mob-suite"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "Mob-suite documentation: https://github.com/phac-nml/mob-suite"
    echo
}

## Option defaults
genome_fa=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
    case "${flag}" in
        i) genome_fa="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------
## Check input
[[ ! -f "$genome_fa" ]] && echo "## ERROR: Input file (-i) $genome_fa does not exist" >&2 && exit 1

## Software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/mob-suite-3.1.0

## Bash strict mode
set -euo pipefail

## Report
echo "## Starting script mob-suite.sh"
date
echo
echo
echo "## Genome FASTA file:                    $genome_fa"
echo "## Output dir:                           $outdir"
echo -e "--------------------\n"


# RUN -------------------------------------------------------------------------
## Make output dir
mkdir -p "$outdir"/logs

echo "## Now running Mob-suite..."
mob_recon \
    --infile "$genome_fa" \
    --outdir "$outdir" \
    --num_threads "$SLURM_CPUS_PER_TASK" \
    --force \
    $more_args

## To get the databases, run:
# mob_init
# DB in /fs/ess/PAS0471/jelmer/conda/mob-suite-3.1.0/lib/python3.8/site-packages/mob_suite/databases


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script mob-suite.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
