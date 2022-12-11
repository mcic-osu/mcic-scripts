#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=sourmash_compare
#SBATCH --output=slurm-sourmash-compare-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "           RUN SOURMASH COMPARE (OPTIONALLY WITH ANI)"
    echo "======================================================================"
    echo
    echo "DESCRIPTION:"
    echo "  This script will:"
    echo "    (1) Create sourmash signatures from FASTA files"
    echo "    (2) Rename the signature to get rid of the extension for plotting"
    echo "    (3) Compare the signature with 'sourmash compare'"
    echo "    (4) Plot a dendrogram and distance/similarity matrix with 'sourmash plot'"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <dir>   Input dir with FASTA files"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --ani                   Use ANI (Average Nucleotide Identity) as the similarity metric"
    echo "  --kmer_size     <int>   Kmer size                                   [default: 31]"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Sourmash compare"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Sourmash Compare and exit"
    echo "  -v/--version            Print the version of Sourmash Compare and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/refgenomes -o results/sourmash -k 29"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://sourmash.readthedocs.io"
    echo "  - Paper:  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/"
    echo
}


## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/project/PAS0471/jelmer/conda/sourmash
    set -u
}

## Print version
Print_version() {
    set +e
    Load_software
    sourmash --version
    set -e
}

## Print help for the focal program
Print_help_program() {
    Load_software
    sourmash compare --help
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

## Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Account (project):    $SLURM_JOB_ACCOUNT"
    echo "Job ID:               $SLURM_JOB_ID"
    echo "Job name:             $SLURM_JOB_NAME"
    echo "Memory (per node):    $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

## Set the number of threads/CPUs
Set_threads() {
    set +u
    if [[ "$slurm" = true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            threads="$SLURM_NTASKS"
        else 
            echo "WARNING: Can't detect nr of threads, setting to 1"
            threads=1
        fi
    else
        threads=1
    fi
    set -u
}

## Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
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
kmer_size=31
ani=false && ani_arg=""
prefix="compare"           # Output filename prefix
indir=""
outdir=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --kmer-size )       shift && kmer_size=$1 ;;
        --ani )             ani=true && ani_arg="--ani" ;;
        --more-args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

## Bash script settings
set -euo pipefail

## Determine output file names
[[ "$ani" = true ]] && prefix=ani
csv_out="$outdir"/output/"$prefix".csv     # distance/ANI matrix in CSV format
cmp_out="$outdir"/output/"$prefix".cmp     # distance/ANI matrix in Python format for sourmash plotting

## Create array with input FASTA files
mapfile -t fastas < <(find "$indir" -iname '*fasta' -or -iname '*fa' -or -iname '*fna' -or -iname '*fna.gz')

## Check input
[[ "$indir" = "" ]] && Die "Please specify an input dir with -i/--indir" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"

## Report
echo
echo "=========================================================================="
echo "                STARTING SCRIPT SOURMASH_COMPARE.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Dir with input FASTA files:       $indir"
echo "Output dir:                       $outdir"
echo "Kmer size:                        $kmer_size"
echo "Run ANI analysis:                 $ani"
[[ $more_args != "" ]] && echo "Other arguments for Sourmash compare:$more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Number of FASTA files:            ${#fastas[@]}"
echo "Listing the input FASTA file(s):"
for fasta in "${fastas[@]}"; do ls -lh "$fasta"; done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create output dirs
mkdir -p "$outdir"/signatures "$outdir"/sig_renamed "$outdir"/output "$outdir"/logs

## Create a signature for each FASTA file
echo "# Creating sourmash signatures..."
${e}sourmash sketch dna \
    -p k="$kmer_size" \
    --outdir "$outdir"/signatures \
    "${fastas[@]}"

## Rename signatures so they have short names
echo -e "--------------------\n"
echo "# Renaming sourmash signatures..."
for sig in "$outdir"/signatures/*sig; do
    newname=$(basename "$sig" .fasta.sig | sed 's/^Spades//')
    newfile="$outdir"/sig_renamed/"$(basename "$sig")"

    ${e}sourmash signature rename \
        "$sig" \
        "$newname" \
        -o "$newfile"
done

## Compare genomes
echo -e "--------------------\n"
echo "## Comparing genomes..."
${e}Time sourmash compare \
    -k"$kmer_size" \
    $ani_arg \
    --from-file <(ls "$outdir"/sig_renamed/*) \
    --csv "$csv_out" \
    $more_args \
    -o "$cmp_out"

## Make plots - run separately for PNG and PDF output
echo -e "\n# Creating plots..."
${e}sourmash plot --labels --output-dir "$outdir"/output "$cmp_out"
${e}sourmash plot --labels --pdf --output-dir "$outdir"/output "$cmp_out"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/output/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
