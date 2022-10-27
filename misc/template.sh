#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=canu
#SBATCH --output=slurm-canu-%j.out


# FUNCTIONS --------------------------------------------------------------------
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "            $0: TODO FUNCTION OF THIS SCRIPT"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i FILE           Input"
    echo "    -o DIR            Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -a STRING         Other argument(s) to pass to TODO_THIS_SOFTWARE"
    echo
    echo "UTILITY OPTIONS"
    echo "    -h                Print this help message and exit"
    echo "    -v                Print the version of TODO_THIS_SOFTWARE and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "    sbatch $0 ..."
    echo
    echo "HARDCODED PARAMETERS:"
    echo "    - ..."
    echo
    echo "OUTPUT:"
    echo "    - ..."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - ..."
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    source activate TODO_THIS_SOFTWARE_ENV
}

## Print version
Print_version() {
    Load_software
    TODO_THIS_SOFTWARE --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}

# CONSTANT AND OPTION DEFAULTS -------------------------------------------------
## Constants

## Option defaults

## Placeholder defaults
infile=""
outdir=""
#tree="" && tree_arg=""

## Parse command-line options
while getopts ':i:o:a:vh' flag; do
    case "${flag}" in
        i) infiles="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        v) Print_version; exit 0 ;;
        h) Print_help; exit 0 ;;
        \?) Print_help; Die "Invalid option $OPTARG" ;;
        :) Print_help; Die "Option -$OPTARG requires an argument" ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
Load_software

## Bash script settings
set -euo pipefail

## Check input
[[ $infile = "" ]] && Die "Please specify an input file with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"


## Report
echo "=========================================================================="
echo "               STARTING SCRIPT TODO_SCRIPTNAME"
date
echo "=========================================================================="
echo "Input file:                  $infiles"
echo "Output dir:                           $outdir"
[[ $more_args != "" ]] && echo "Other arguments for TODO_THIS_SOFTWARE:    $more_args"
echo "Listing input file:"
ls -lh "$infile"
echo "=========================================================================="


# MAIN -------------------------------------------------------------------------
## Create the output directory
mkdir -p "$outdir"/logs

echo -e "\n## Now running TODO_THIS_SOFTWARE..."
set -o xtrace
TODO_COMMAND \
    -t "$SLURM_CPUS_PER_TASK"
    $more_args \
set +o xtrace


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
