#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name=gubbins
#SBATCH --output=slurm-gubbins-%j.out

#TODO - Process dates once Gubbins 3.3 comes on conda -- see `--date` argument in docs

# FUNCTIONS --------------------------------------------------------------------
## Help function
print_help() {
    echo
    echo "======================================================================"
    echo "$0: Run Gubbins to identify recombination among bacterial genomes"
    echo "======================================================================"
    echo "USAGE:"
    echo "------------------"
    echo "$0 -i <alignment-file> -o <output-prefix> ..."
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i FILE           Input alignment in FASTA format"
    echo "    -o STRING         Output prefix (dir + run ID)"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -t FILE           Tree file for starting tree"
    echo "    -a STRING         Other argument(s) to pass to Gubbins"
    echo "    -v                Print the version and exit"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "sbatch $0 -i aln.fasta -o results/gubbins/run1"
    echo
    echo "DOCUMENTATION:"
    echo "------------------"
    echo "Documentation:        http://nickjcroucher.github.io/gubbins/ "
    echo "Tutorial:             https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_tutorial.md "
    echo "Manual:               https://github.com/nickjcroucher/gubbins/blob/master/docs/gubbins_manual.md "
    echo "Repository:           https://github.com/nickjcroucher/gubbins "
    echo "Paper:                https://academic.oup.com/nar/article/43/3/e15/2410982 "
    echo
}

## Load software
load_software() {
    module load python/3.6-conda5.2
    source activate /fs/project/PAS0471/jelmer/conda/gubbins-3.2.1
}

## Print the version
print_version() {
    load_software
    run_gubbins.py --version
}


# PARSE OPTIONS ----------------------------------------------------------------
## Option defaults
aln=""
tree="" && tree_arg=""
out_prefix=""
more_args=""

## Parse command-line options
while getopts ':i:t:o:a:vh' flag; do
    case "${flag}" in
        i) aln="$OPTARG" ;;
        o) out_prefix="$OPTARG" ;;
        t) tree="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        v) print_version && exit 0 ;;
        h) print_help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ "$aln" = "" ]] && echo "## ERROR: Please specify an alignment with -i" >&2 && exit 1
[[ "$out_prefix" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$aln" ]] && echo "## ERROR: Input file $aln does not exist" >&2 && exit 1
[[ "$tree" != "" && ! -f "$tree" ]] && echo "## ERROR: Input file $tree does not exist" >&2 && exit 1

## Build tree and dates arguments
[[ $tree != "" ]] && tree_arg="--starting-tree $tree"

## Report
echo "=========================================================================="
echo "                       STARTING SCRIPT GUBBINS.SH"
date
echo "=========================================================================="
echo "## Alignment input file:                      $aln"
echo "## Output prefix:                             $out_prefix"
[[ "$tree" != "" ]] && echo "## Tree input file:                           $tree"
[[ $more_args != "" ]] && echo "## Other arguments for Gubbins:              $more_args"
echo
echo "## Listing the input files:"
ls -lh "$aln"
[[ $tree != "" ]] && ls -lh "$tree"
echo "=========================================================================="


# MAIN -------------------------------------------------------------------------
## Create output dir
outdir=$(dirname "$out_prefix")
mkdir -p "$outdir"/logs

## Run treetime
echo "## Now running Gubbins..."

set -o xtrace
run_gubbins.py \
    --threads "$SLURM_CPUS_PER_TASK" \
    --prefix "$out_prefix" \
    $tree_arg \
    "$aln"
set +o xtrace


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Gubbins version used:"
run_gubbins.py --version
echo
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo 
echo "## Done with script gubbins.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
