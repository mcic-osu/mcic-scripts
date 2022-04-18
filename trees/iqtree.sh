#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-iqtree-%j.out


# ARGS AND PARAMETERS ----------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Construct a tree from a FASTA alignment"
    echo
    echo "Syntax: $0 -i <input-FASTA> -o <output dir> ..."
    echo
    echo "Required options:"
    echo "  -i FILE       Input FASTA file -- should contain multiple, aligned sequences"
    echo "  -p STRING     Output prefix: directory (will be created if needed) + first part of filename"
    echo
    echo "Other options:"
    echo "  -a STRING     Other arguments to pass to IQ-tree"
    echo "  -b NUM        Number of bootstraps                  [default: no bootstrapping]"
    echo "  -h            Print this help message and exit"
    echo
    echo "Example command:"
    echo "$0 -i results/alignment/COI_aligned.fa -p results/iqtree/COI -b 1000"
    echo
    echo "To submit to the OSC queue, preface with sbatch:"
    echo "sbatch $0 ..."
    echo
    echo "SLURM parameters in script: '--account=PAS0471 --time=3:00:00 --cpus-per-task=12 --mem=32G --output=slurm-iqtree-%j.out"
    echo
    echo "Default SLURM parameters can be overridden when submitting the script, e.g.:"
    echo "sbatch -t 15 $0 ...      (override default time reservation of 3 hours, use 15 minutes instead)"
    echo
    echo "IQ-tree documentation: http://www.iqtree.org/doc/"
    echo
}

## Option defaults
fa_in=""
prefix=""
more_args=""
nboot=""

## Parse options
while getopts ':i:p:b:a:h' flag; do
    case "${flag}" in
        i) fa_in="$OPTARG" ;;
        p) prefix="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        b) nboot="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo "## ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
        :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Report
echo "## Starting script iqtree.sh"
date
echo

## Check input
[[ "$fa_in"  = "" ]] && echo "## ERROR: Please provide input FASTA file with -i flag" && exit 1
[[ "$prefix"  = "" ]] && echo "## ERROR: Please provide output prefix with -p flag" && exit 1
[[ ! -f $fa_in ]] && echo "## ERROR: Input FASTA file ($fa_in) does not exist" && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/iqtree-2.2.0

## Bash strict settings
set -euo pipefail

## Other parameters
n_cores="$SLURM_CPUS_ON_NODE"

if [[ "$nboot" != "" ]]; then
    boot_arg="-B $nboot"
else
    boot_arg=""
fi

## Report
echo "## Input FASTA:                  $fa_in"
echo "## Output prefix:                   $prefix"
[[ "$nboot" != "" ]] && echo "## Number of bootstraps:                   $nboot"
[[ "$more_args" != "" ]] && echo "## Other arguments to pass to IQ-tree:                   $more_args"
echo -e "-----------------------------\n"

## Create output dir
mkdir -p "$(dirname "$prefix")"


# RUN IQ-TREE --------------------------------------------------------------
echo -e "## Starting IQ-tree run..."
iqtree \
    -s "$fa_in" \
    --prefix "$prefix" \
    -redo \
    -m MFP \
    -nt AUTO -ntmax "$n_cores" $boot_arg $more_args

#? -m MFP  => Model selection with ModelFinder (is the IQ-tree default, too)
#? -redo   => Will overwrite old results
#? --mem 4G  - memory

# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$prefix"*
echo -e "\n## Done with script iqtree.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
