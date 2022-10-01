#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=krone
#SBATCH --output=slurm-krona-%j.out

# HELP AND COMMAND-LINE OPTIONS ------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Make Krona plots of Kraken output"
    echo
    echo "Syntax: $0 -i <input-file> -o <output-file>"
    echo 
    echo "Required options:"
    echo "  -i FILE           Input file (a Kraken2 '_main.txt' output file)"
    echo "  -o FILE           Output file (will be in HTML format)"
    echo
    echo "Other options:"
    echo "  -h                Print this help message and exit"
    echo
    echo "Example:          $0 -i results/kraken/sampleA_main.txt -o results/krona/sampleA.html"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
infile=""
outfile=""

## Get command-line options
while getopts 'i:o:h' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outfile="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/krona-env

## Bash strict settings
set -euo pipefail

## Check input options
[[ "$infile" = "" ]] && echo "ERROR: Must specify input file with -i" >&2 && exit 1
[[ "$outfile" = "" ]] && echo "ERROR: Must specify input file with -i" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "ERROR: Input file $infile does not exist" >&2 && exit 1

## Make outdir only if outfile is not in current dir (contains a "/")
if echo "$outfile" | grep -q "/"; then
    outdir=$(dirname "$outfile")
    mkdir -p "$outdir"
fi

## Report
echo "## Starting script krona.sh..."
date
echo
echo "## Input file:           $infile"
echo "## Output file:          $outfile"
echo -e "----------------------------\n"


# MAIN -------------------------------------------------------------------------
ktImportTaxonomy -q 2 -t 3 "$infile" -o "$outfile"

# q: column with query ID
# t: column with taxonomy ID


# WRAP UP ----------------------------------------------------------------------
echo
echo "## Listing output file:"
ls -lh "$outfile"
echo
echo "## Done with script krona.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
