#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=concat_assemblies
#SBATCH --output=slurm-concat_assemblies-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              CONCATENATE TRANSCRIPTOME ASSEMBLIES"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -o <output file> <assembly 1> <assembly 2> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile     <dir>   Output dir (will be created if needed)"
    echo "  One or more assembly FASTA files as positional arguments"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -o results/concat_assemblies.fasta results/spades.fasta results/trinity.fasta"
    echo
    echo "OUTPUT:"
    echo "  - Each contig will have an added ID from the filename of the original assembly"
    echo
}

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo >&2
    echo "=====================================================================" >&2
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option" >&2
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'" >&2
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    echo -e "\nEXITING..." >&2
    echo "=====================================================================" >&2
    echo >&2
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
declare -a assemblies
outdir=""

## Parse command-line args
all_args="$*"
count=0

while [ "$1" != "" ]; do
    case "$1" in
        -o | --outfile )    shift && outfile=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h | --help )       Print_help; exit 0 ;;
        --dryrun )          dryrun=true ;;
        --debug )           debug=true ;;
        * )                 assemblies[$count]=$1 && count=$(( count + 1 )) ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Bash script settings
set -euo pipefail

## Find the output dir
outdir=$(dirname "$outfile")

## Check input
[[ ${#assemblies[@]} = 0 ]] && Die "Please specify one or more input assemblies" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT CONCAT_ASSEMBLIES.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Output file:                      $outfile"
echo
echo "Listing the input assemblies:"
for assembly in "${assemblies[@]}"; do
    [[ ! -f $assembly ]] && Die "Input assembly $assembly does not exist!"
    ls -lh "$assembly"
done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    ## Create the output dir
    mkdir -p "$outdir"/logs

    true > "$outfile"

    for assembly in "${assemblies[@]}"; do
        asm_id=$(basename "$assembly" .fasta)
        sed "s/>/>${asm_id}_/" "$assembly" >> "$outfile"
    done

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outfile"
fi
echo "# Done with script"
date
