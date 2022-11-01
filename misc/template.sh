#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=TODO_THIS_SOFTWARE
#SBATCH --output=slurm-TODO_THIS_SOFTWARE-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
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
    echo "    -x                Run the script in debug mode"
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

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Constants

## Option defaults
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
indir=""
outdir=""
#tree="" && tree_arg=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        -X | --debug )          debug=true ;;
        -N | --dryrun )         dryrun=false ;;
        -v | --version )        Print_version && exit ;;
        -h | --help )           Print_help && exit ;;
        * )                     Die "Invalid option $1" && exit 1 ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## Bash script settings
set -euo pipefail

## Check input
[[ $indir = "" ]] && Die "Please specify an input dir with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -d $indir ]] && Die "Input dir $indir does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT TODO_SCRIPTNAME"
date
echo "=========================================================================="
echo "Input dir:                   $indir"
echo "Output dir:                  $outdir"
[[ $more_args != "" ]] && echo "Other arguments for TODO_THIS_SOFTWARE:    $more_args"
echo "Listing input file:"
ls -lh "$indir"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs
fi

echo -e "\n## Now running TODO_THIS_SOFTWARE..."
[[ "$debug" = false ]] && set -o xtrace

TODO_COMMAND \
    -t "$SLURM_CPUS_PER_TASK"
    $more_args \

[[ "$debug" = false ]] && set +o xtrace


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "## Done with script"
date
if [[ "$dryrun" = false ]]; then
    echo -e "\n## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
