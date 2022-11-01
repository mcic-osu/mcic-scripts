#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --job-name=fqconcat
#SBATCH --out=slurm-fqconcat-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "      $0: Concatenate FASTQ files from different lanes"
    echo "======================================================================"
    echo
    echo "PURPOSE:"
    echo "  This script will concatenate FASTQ files from different lanes"
    echo "  (but the same sample and read direction)."
    echo "  NOTE: The script assumes that files for each sample are in their own subdirectory." 
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i DIR         Input dir containing FASTQ files, in separate underlying dirs (1 per sample)"
    echo "  -o DIR         Output dir for concatenated FASTQ files (all will be placed directly in this dir)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -n INTEGER     Number of files expected per sample per read direction       [default: 2]"
    echo "  -p STRING      Globbing pattern for underlying dirs                         [default: '*' i.e. any dir]"
    echo "  -r STRING      Literal string or regex pattern to remove from file names    [default: none]"
    echo
    echo "UTILITY OPTIONS"
    echo "  -h             Print this help message and exit"
    echo "  -d             Dryrun (don't execute commands)"
    echo "  -x             Run the script in debug mode (print all code)"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/original -o data/fastq/concat"
    echo "  sbatch $0 -i data/fastq/original -o data/fastq/concat -n 2"
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
nfiles=2            # Expect two files per lane
se=false            # Single-end=false, i.e. expect paired-end seqs
debug=false
dryrun=false

indir=""
outdir=""
dir_pattern="*"
remove_pattern=""

## Parse command-line options
while getopts ':i:o:p:r:n:sdhx' flag; do
    case "${flag}" in
        i) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        p) dir_pattern="$OPTARG" ;;
        r) remove_pattern="$OPTARG" ;;
        n) nfiles="$OPTARG" ;;
        s) se=true ;;
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

## Input checks
[[ $indir = "" ]] && Die "Please specify an input dir with -i"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ "$indir" = "$outdir" ]] && Die "Input dir $indir should not be the same as the output dir"
[[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FQCONCAT.SH"
date
echo "=========================================================================="
echo "Input dir:                                  $indir"
echo "Output dir:                                 $outdir"
[[ "$dir_pattern" != "" ]] && echo "Dir globbing pattern:                       $dir_pattern"
[[ "$remove_pattern" != "" ]] && echo "Pattern/string to remove from file names:   $remove_pattern"
echo "Expected nr files for each sample & read direction:   $nfiles"
[[ "$dryrun" = true ]] && echo "THIS IS A DRY RUN"
echo "=========================================================================="
echo

# ==============================================================================
#                               RUN
# ==============================================================================
## Create output dir, if necessary
[[ "$dryrun" = false ]] && mkdir -p "$outdir"

## Find per-sample directories
dirs=( $(find $indir -mindepth 1 -type d -name "*$dir_pattern*") )
echo "## Nr directories (samples) found:          ${#dirs[@]}"
[[ ${#dirs[@]} -eq 0 ]] && echo "Die: No matching dirs found!"
echo

## Loop over per-sample directories
for fdir in "${dirs[@]}"; do
    
    ## Count input files - exit if there are not 2 files per read direction
    n_R1=$(find "$fdir" -maxdepth 1 -name "*_R1[._]*fastq.gz" | wc -l)
    [[ "$se" = false ]] && n_R2=$(find "$fdir" -maxdepth 1 -name "*_R2[._]*fastq.gz" | wc -l)

    ## Ignore dirs with no FASTQ files -- e.g. higher level dirs (!)
    if [ "$n_R1" != 0 ]; then
        echo -e "\n## Directory: $fdir"

        ## Check input
        [[ "$n_R1" != "$nfiles" ]] && Die "Did not find $nfiles R1 input files"
        [[ "$se" = false ]] && [[ "$n_R2" != "$nfiles" ]] && Die "Did not find $nfiles R2 input files"

        ## Each dir should contain FASTQ files for a single sample.
        ## Extract the sample ID (e.g. "S9") from the 1st FASTQ file name.
        smp_id=$(basename "$(find "$fdir" -maxdepth 1 -name "*fastq.gz" | head -n 1)" .fastq.gz)
        smp_id=$(echo "$smp_id" | sed -E 's/_L0[0-9][0-9]_.*//')
        smp_id="${smp_id/$remove_pattern/}"
        echo "## Sample ID: $smp_id"

        ## Make input dir read-only
        [[ "$dryrun" = false ]] && chmod a-w "$fdir"/*fastq.gz

        ## List input files
        echo "## R1 input files:"
        find "$fdir" -name "*_R1[._]*.fastq.gz" -print0 | sort -z | xargs -0 du -sh
        [[ "$se" = false ]] && echo "## R2 input files:"
        [[ "$se" = false ]] && find "$fdir" -name "*_R2[._]*.fastq.gz" -print0 | sort -z | xargs -0 du -sh
        
        ## Output files
        R1_out="$outdir"/"$smp_id"_R1_001.fastq.gz
        [[ "$se" = false ]] && R2_out="$outdir"/"$smp_id"_R2_001.fastq.gz

        if [[ "$dryrun" = false ]]; then
            ## Concatenate FASTQ files
            find "$fdir" -name "*_R1[._]*.fastq.gz" -print0 | sort -z | \
                xargs -0 -I{} cat {} > "$R1_out"
            
            [[ "$se" = false ]] && find "$fdir" -name "*_R2[._]*.fastq.gz" -print0 | \
                sort -z | xargs -0 -I{} cat {} > "$R2_out"

            ## List output files
            echo "## Output files:"
            du -h "$R1_out"
            [[ "$se" = false ]] && du -h "$R2_out"
        else
            echo "## Output files:"
            echo "$R1_out"
            [[ "$se" = false ]] && echo "$R2_out"
        fi
    fi
    
    echo -e "------------------\n"
done

## Make concatenated files read-only to protect against accidental removal and overwriting
echo "## Making concatenated files read-only..."
chmod a-w "$outdir"/*fastq.gz


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
