#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=nanoplot
#SBATCH --output=slurm-nanoplot-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "      $0: Run Nanoplot to QC ONT/PACBIO READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--infile FILE      Input file: sequencing summary from Albacore or Guppy"
    echo "    -o/--outdir DIR       Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -a/--more_args \"QUOTED STRING\"   Other argument(s) to pass to Nanoplot"
    echo
    echo "UTILITY OPTIONS"
    echo "    -x/--debug            Run the script in debug mode"
    echo "    -h/--help             Print this help message and exit"
    echo "    -v/--version          Print the version of Nanoplot and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "    sbatch $0 -i results/guppy/sequencing_summary.txt -o results/nanoplot"
    echo "    sbatch $0 -i results/guppy/sequencing_summary.txt -o results/nanoplot -a \"--minlength N\"" 
    echo
    echo "HARDCODED PARAMETERS:"
    echo "    - This script assumes that the input file is a sequencing summary file,"
    echo "      whereas Nanoplot can also be run on actual sequence files."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - https://github.com/wdecoster/NanoPlot"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/nanoplot
}

## Print version
Print_version() {
    Load_software
    NanoPlot --version
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
        -N | --dryrun )         dryrun=false ;;
        -v | --version )        Print_version && exit ;;
        -h | --help )           Print_help && exit ;;
        * )                     Die "Invalid option $1" ;;
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
[[ $infile = "" ]] && Die "Please specify an input file with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input dir $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT NANOPLOT.SH"
date
echo "=========================================================================="
echo "Input file:                  $infile"
echo "Output dir:                  $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Nanoplot:    $more_args"
echo "Listing input file:"
ls -lh "$infile"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs
fi

echo -e "\n## Now running Nanoplot..."
[[ "$debug" = false ]] && set -o xtrace

NanoPlot \
    --summary "$infile" \
    --outdir "$outdir" \
    --tsv_stats \
    --info_in_report \
    --threads "$SLURM_CPUS_PER_TASK" \
    $more_args

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
    echo
fi
echo "## Done with script"
date
