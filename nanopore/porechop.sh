#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --job-name=porechop
#SBATCH --output=slurm-porechop-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "             RUN PORECHOP TO REMOVE ONT ADAPTERS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <dir>   Input dir"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -a/--more_args  <string> Quoted string with additional argument(s) to pass to Porechop"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -N/--dryrun             Dry run: don't execute commands, only parse arguments and report"
    echo "  -x/--debug              Run the script in debug mode (print all code)"
    echo "  -v/--version            Print the version of Porechop and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i TODO -o results/TODO "
    echo "  sbatch $0 -i TODO -o results/TODO -a \"-x TODO\""
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
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate; done
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/porechop-0.2.4
}

## Print version
Print_version() {
    Load_software
    porechop --version
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
## Option defaults
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
outdir=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        -X | --debug )          debug=true ;;
        -N | --dryrun )         dryrun=true ;;
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
if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
    threads="$SLURM_CPUS_PER_TASK"
elif [[ -n "$SLURM_NTASKS" ]]; then
    threads="$SLURM_NTASKS"
else
    threads=1
fi

## Bash script settings
set -euo pipefail

## Output file
outfile="$outdir"/$(basename "$infile")

## Check input
[[ $infile = "" ]] && Die "Please specify an input file with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT PORECHOP.SH"
date
echo "=========================================================================="
echo "Input file:                  $infile"
echo "Output dir:                  $outdir"
echo
echo "Output file                  $outfile"
[[ $more_args != "" ]] && echo "Other arguments for Porechop:    $more_args"
echo "Listing input file:"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN\n"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    
    ## Create the output directory
    mkdir -p "$outdir"/logs

    ## Run
    echo -e "\n## Now running Porechop..."
    [[ "$debug" = false ]] && set -o xtrace
    porechop \
        -i "$infile" \
        -o "$outfile" \
        --threads "$threads"
        $more_args \
    [[ "$debug" = false ]] && set +o xtrace

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    echo
    echo "========================================================================="
    echo "## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
