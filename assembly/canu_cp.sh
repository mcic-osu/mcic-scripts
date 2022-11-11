#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=canu_cp
#SBATCH --output=slurm-canu_cp-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Print_help() {
    echo
    echo "==============================================================================="
    echo "                $0: Copy Cany output"
    echo "==============================================================================="
    echo "Because the main Canu job stops well before all its SLURM jobs are done,"
    echo "there is now way to copy its output files in the main script."
    echo "Additionally, because Canu produces so much output,"
    echo "we want the initial output dir to be a scratch dir and then copy selected files." 
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <canu workdir> -o <assembly outfile>"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  (Choose one of the input file options: -i or -I)"
    echo "  -i/--indir      <dir>  The scratch dir with all Canu output"
    echo "  -o/--outfile    <file> The final output dir for selected Canu output files"
    echo
    echo "UTILITY OPTIONS"
    echo "  --dryrun               Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                Run the script in debug mode (print all code)"
    echo "  -h                     Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i /fs/scratch/PAS0471/canu -o results/canu/assembly.fasta"
    echo
}

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Defaults
dryrun=false
debug=false

## Placeholder defaults
indir=""
outfile=""

## Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -o | --outfile )        shift && outfile=$1 ;;
        -h )                    Print_help; exit 0 ;;
        --dryrun )              dryrun=true && e="echo " ;;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Bash strict settings
#set -euo pipefail

## Check input
[[ $indir = "" ]] && Die "Please specify input dir with -i or -i" "$all_args"
[[ $outfile = "" ]] && Die "Please specify an output file with -o" "$all_args"

## Determine output dir
outdir=$(dirname "$outfile")

## Report
echo "=========================================================================="
echo "                      STARTING SCRIPT CANU_CP.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo
echo "Input dir:                        $indir"
echo "Output asssembly file:            $outfile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    ## Create the output directory
    mkdir -p "$outdir"/logs

    ## Copy files to final outdir
    echo -e "\n# Copying files to final outdir..."
    cp -v "$indir"/* "$outdir"    # Will only copy files, not dirs - those are huge
    cp -r "$indir"/canu-logs "$outdir"

    ## Copy assembly FASTA
    echo -e "\n# Copying the assembly FASTA file:"
    cp -v "$indir"/*.contigs.fasta "$outfile"
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing the final assembly file:"
    ls -lh "$outfile"
fi
echo -e "\n# Done with script"
date
