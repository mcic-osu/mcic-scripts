#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --job-name=fqconcat
#SBATCH --out=slurm-fqconcat-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                      $0"
    echo "              Concatenate FASTQ files from different lanes"
    echo "======================================================================"
    echo
    echo "PURPOSE:"
    echo "  This script will concatenate gzipped FASTQ files"
    echo "    - Usually, this will involve concatenating files from the sample sample but different lanes."
    echo "    - Files _can_ be across different directories, but if you have one dir per sample,"
    echo "      it might be easier to run the 'fq_concat.sh' script where you don't need a samplesheet."
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <samplesheet> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i FILE        Samplesheet with two (tab-separated) columns (tab-separated, no header):"
    echo "                   - Column 1 should have filenames and column 2 should have IDs (sample IDs/groupings)"
    echo "                   - All files with the same ID in the second column will be concatenated"
    echo "                   - For paired-end files, include the read direction in the ID"
    echo "                   - The names of output files will be: <outdir>/<sampleID>.fastq.gz"
    echo "  -o DIR         Output dir for concatenated FASTQ files"
    echo "                   - All will be placed directly in this dir"
    echo "                   - This should not be the same dir as (one of) the input dir(s)"
    echo
    echo "UTILITY OPTIONS"
    echo "  -h             Print this help message and exit"
    echo "  -d             Dryrun (don't execute commands)"
    echo "  -x             Run the script in debug mode (print all code)"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/meta/samples.txt -o data/fastq/concat"
    echo
    echo "EXAMPLE SAMPLESHEET FOR PAIRED-END FILES:"
    echo "data/fastq/sampleA_L001_R1.fastq.gz   sampleA_R1"
    echo "data/fastq/sampleA_L001_R2.fastq.gz   sampleA_R2"
    echo "data/fastq/sampleA_L002_R1.fastq.gz   sampleA_R1"
    echo "data/fastq/sampleA_L002_R2.fastq.gz   sampleA_R2"
    echo "data/fastq/sampleB_L001_R1.fastq.gz   sampleB_R1"
    echo "data/fastq/sampleB_L001_R2.fastq.gz   sampleB_R2"
    echo "data/fastq/sampleB_L002_R1.fastq.gz   sampleB_R1"
    echo "data/fastq/sampleB_L002_R2.fastq.gz   sampleB_R2"
    echo
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Defaults
debug=false
dryrun=false

samplesheet=""
outdir=""

## Parse command-line options
while getopts ':i:o:dxh' flag; do
    case "${flag}" in
        i) samplesheet="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) dryrun=true ;;
        x) debug=true ;;
        h) Print_help && exit 0 ;;
        \?) Die "Invalid option -$OPTARG" ;;
        :) Die "Option -$OPTARG requires an argument" ;;
    esac
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Bash strict settings
set -euo pipefail

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FQCONCAT.SH"
date
echo "=========================================================================="
echo "Sample sheet:             $samplesheet"
echo "Output dir:               $outdir"
[[ "$dryrun" = true ]] && echo "THIS IS A DRY RUN"
echo "=========================================================================="
echo

## Input checks
[[ $samplesheet = "" ]] && Die "Please specify a sample sheet with -i"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f "$samplesheet" ]] && Die "Sample sheet $samplesheet does not exist"


# ==============================================================================
#                               RUN
# ==============================================================================
## Create output dir, if necessary
[[ "$dryrun" = false ]] && mkdir -p "$outdir"

## Loop over IDs, concatenate:
for ID in $(cut -f 2 "$samplesheet" | sort | uniq); do

    outfile="$outdir"/"$ID".fastq.gz
    infiles=( "$(awk -v smp="$ID" '$2 == smp' "$samplesheet" | cut -f 1)" )

    echo "## Input files for sample $ID:"
    ls -lh $infiles
    
    ## Check that indir isn't the same as the outdir
    indir=$(dirname $infiles | head -n 1)
    [[ "$indir" = "$outdir" ]] && "Die: Input dir should not be the same as the output dir"

    [[ "$dryrun" = false ]] && cat $infiles > "$outfile"

    echo "## Output file:"
    [[ "$dryrun" = false ]] && du -h "$outfile"
    [[ "$dryrun" = true ]] && echo "$outfile"
    echo -e "-----------------\n"
done

## Make concatenated files read-only to protect against accidental removal and overwriting
if [[ "$dryrun" = false ]]; then
    echo "## Making concatenated files read-only..."
    chmod a-w "$outdir"/*fastq.gz
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo -e "\n====================================================================="
echo -e "## Done with script"
date
if [[ "$dryrun" = false ]]; then
    echo "## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
