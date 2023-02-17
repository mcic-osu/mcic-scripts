#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=deeptmhmm
#SBATCH --output=slurm-deeptmhmm-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run DeepTMHMM."
    echo
    echo "Syntax: $0 -i <input-FASTA> -o <output-dir>..."
    echo
    echo "Required options:"
    echo      "-i STRING         Input amino acid FASTA file"
    echo      "-o STRING         Output directory"
    echo
    echo "Other options:"
    echo      "-h                Print this help message and exit"
    echo
    echo "Example:               $0 -i ref/aa.fa -o results/tmhmm"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "DeepTMHMM paper: https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1"
    echo "DeepTMHMM online server: https://dtu.biolib.com/DeepTMHMM"
    echo
}

## Option defaults
infile=""
outdir=""

## Parse command-line options
while getopts ':i:o:h' flag; do
    case "${flag}" in
        i) infile="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done

## Check input
[[ "$infile" = "" ]] && echo "## ERROR: Please specify an input FASTA file with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output directory with -o" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input FASTA file (-i) $infile does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/deeptmhmm2
DEEPTMHMM_DIR=/fs/ess/PAS0471/jelmer/conda/deeptmhmm2/biolib_container/openprotein/

## Bash strict mode
set -euo pipefail

## Make input and output absolute
[[ ! "$infile" =~ ^/ ]] && infile="$PWD"/"$infile"
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"

## Report
echo
echo "Starting script deeptmhmm.sh"
date
echo "Input AA FASTA:    $infile"
echo "Output dir:        $outdir"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
## Make output dir
mkdir -p "$outdir"

echo "## Now running DeepTMHMM..."

cd "$DEEPTMHMM_DIR" || exit

python3 predict.py \
    --fasta "$infile" \
    --outdir "$outdir"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script deeptmhmm.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
