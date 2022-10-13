#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=16
#SBATCH --job-name=roary
#SBATCH --time=6:00:00
#SBATCH --output=slurm-roary-%j.output

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run Roary for a pangenome analysis"
    echo
    echo "Syntax: $0 -i <input-dir> -o <output-dir>..."
    echo
    echo "Required options:"
    echo "    -i DIR            Input dir with GFFs produced by Prokka"
    echo "                      NOTE: The GFF files should contain the nucleotide sequences, too."
    echo "    -o DIR            Output dir"
    echo
    echo "Other options:"
    echo "    -a STRING         Other argument(s) to pass to Roary"
    echo "    -h                Print this help message and exit"
    echo
    echo "Example:              $0 -i results/prokka -o results/roary"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "Roary documentation: https://sanger-pathogens.github.io/Roary/"
    echo "Roary paper: http://bioinformatics.oxfordjournals.org/content/31/22/3691"
    echo
}

## Option defaults
indir=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
    case "${flag}" in
        i) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done

# SETUP ------------------------------------------------------------------------
## Check input
[[ "$indir" = "" ]] && echo "## ERROR: Please specify an input dir with -i" >&2 && exit 1
[[ "$outdir" = "" ]]  && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-i) $indir does not exist" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/roary-3.13

## Bash strict settings
set -euo pipefail

## Report
echo "## Starting script roary.sh"
date
echo
echo "Input dir:          $indir"
echo "Outdir dir:         $outdir"
[[ $more_args != "" ]] && echo "## Other arguments to pass to Roary:    $more_args"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
roary \
    -v "$indir"/*gff \
    -f "$outdir" \
    -e -n \
    -p "$SLURM_CPUS_PER_TASK" $more_args

#? -e Creates a multiFASTA alignment of core genes
#? -n aligns core genes with MAFFT


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script roary.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
