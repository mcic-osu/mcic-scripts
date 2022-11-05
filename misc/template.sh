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
    echo "                            $0"
    echo "                  TODO FUNCTION OF THIS SCRIPT"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--indir FILE        Input"
    echo "    -o/--outdir DIR        Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -a/--more_args STRING  Quoted string with additional argument(s) to pass to TODO_THIS_SOFTWARE"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h/--help              Print this help message and exit"
    echo "    -N/--dryrun            Dry run: don't execute commands, only parse arguments and report"
    echo "    -x/--debug             Run the script in debug mode (print all code)"
    echo "    -v/--version           Print the version of TODO_THIS_SOFTWARE and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "    sbatch $0 -i TODO -o results/TODO "
    echo "    sbatch $0 -i TODO -o results/TODO -a \"-x TODO\""
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
    module load miniconda3/4.12.0-py39
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
more_args=""
#tree="" && tree_arg=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        -X | --debug )          debug=true ;;
        -N | --dryrun )         dryrun=false ;;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_help; Die "Invalid option $1" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## Get number of threads
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## FASTQ filename parsing
extension=$(echo "$R1_in" | sed -E 's/.*(\.fa?s?t?q\.gz$)/\1/')
R1_suffix=$(echo "$R1_in" | sed -E "s/.*(_R?1)_?[[:digit:]]*$extension/\1/")
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1_in" | sed -E "s/${R1_suffix}_?[[:digit:]]*${extension}//")

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
    -t "$n_threads"
    $more_args \

[[ "$debug" = false ]] && set +o xtrace


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
