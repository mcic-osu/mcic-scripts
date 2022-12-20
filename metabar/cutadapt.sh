#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=180
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=cutadapt
#SBATCH --output=slurm-cutadapt-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "================================================================================"
    echo "                            $0"
    echo " Run Cutadapt to remove metabarcoding primers for a single pair of FASTQ files"
    echo "================================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-R1-FASTQ> -o <output-dir> [ -f <fwd-primer> -r <rev-primer | -p <primer-file> ] [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Input R1 FASTQ file (corresponding R2 will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  NOTE: You should also provide either a pair of primer sequences or a file with primers sequences, see below"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --primer-f      <str>   Forward primer sequence (use in combination with -r)"
    echo "  --primer-r      <str>   Reverse primer sequence (use in combination with -f)"
    echo "  --primer-file   <file>  File with primer sequences, one pair per line separated by a space (alternative to using -f and -r)"
    echo "  --keep-untrimmed        Don't discard untrimmed sequences (default: discard sequences with no primers)"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Cutadapt"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Cutadapt and exit"
    echo "  -v/--version            Print the version of Cutadapt and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/sample1_R1.fastq.gz -o results/cutadapt -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  The '--pair-filter=any' option is always used."
    echo 
    echo "NOTES:"
    echo "  - When you have multiple primer pairs, specify a primer file with -p."
    echo "  - The script will compute and use the reverse complements of both primers."

    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://cutadapt.readthedocs.io/en/stable/"
    echo "  - Paper: https://journal.embnet.org/index.php/embnetjournal/article/view/200/0"
    echo
}

## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/cutadapt-4.1
    set -u
}

## Print version
Print_version() {
    set +e
    Load_software
    cutadapt --version
    set -e
}

## Print help for the focal program
Print_help_program() {
    Load_software
    cutadapt --help
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
discard_untrimmed=true && discard_arg="--discard-untrimmed"

debug=false
dryrun=false && e=""
slurm=true

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
R1_in=""
outdir=""
primer_f=""
primer_r=""
primer_file=""
primer_arg=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )             shift && R1_in=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --primer-f )            shift && primer_f=$1 ;;
        --primer-r )            shift && primer_r=$1 ;;
        --primer-file )         shift && primer_file=$1 ;;
        --keep-untrimmed )      discard_untrimmed=false ;;
        --more-args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0;;
        --dryrun )              dryrun=true && e="echo ";;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## Bash script settings
set -euo pipefail

## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

## Determine input dir and R2 file
indir=$(dirname "$R1_in")
file_ext=$(basename "$R1_in" | sed -E 's/.*(.fastq|.fq|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1_in" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1_in" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")
R1_basename=$(basename "$R1_in")
R2_basename=$(basename "$R2_in")

## Check input
[[ "$R1_in" = "" ]] && Die "Please specify an R1 input file with -i/--R1" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$R1_in" ]] && Die "Input FASTQ file $R1_in not found"
[[ ! -f "$R2_in" ]] && Die "Input FASTQ file $R2_in not found"
[[ "$indir" = "$outdir" ]] && Die "Input dir should not be the same as output dir"

## Define cutadapt options
[[ "$discard_untrimmed" = "false" ]] && discard_arg=""

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT CUTADAPT.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input FASTQ file R1:              $R1_in"
echo "Output dir:                       $outdir"
echo "Discard untrimmed (-d):           $discard_untrimmed"
[[ $more_args != "" ]] && echo "Other arguments for Cutadapt:$more_args"
echo "Input FASTQ file R2 (inferred):   $R2_in"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$R1_in" "$R2_in"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               DEFINE PRIMERS
# ==============================================================================
## Get primers
if [ "$primer_file" = "" ]; then
    echo "Using forward and reverse primers provided as arguments..."
    [[ $primer_f = "" ]] && Die "No forward primer (-f) provided"
    [[ $primer_r = "" ]] && Die "No reverse primer (-r) provided"

    primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
    primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

    primer_arg="-a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"

    echo "Forward primer (-f):              $primer_f"
    echo "Reverse primer (-r):              $primer_r"
    echo "Forward primer - rev. comp.:      $primer_f_rc"
    echo "Reverse primer - rev. comp.:      $primer_r_rc"

else
    echo "Using primer file $primer_file to read primers..."
    [[ ! -f "$primer_file" ]] && Die "Primer file $primer_file not found"

    while read -r primer_f primer_r; do

        primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
        primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

        echo
        echo "Forward primer (-f):          $primer_f"
        echo "Reverse primer (-r):          $primer_r"
        echo "Forward primer - rev. comp.:  $primer_f_rc"
        echo "Reverse primer - rev. comp.:  $primer_r_rc"

        primer_arg="$primer_arg -a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"
        primer_arg=$(echo "$primer_arg" | sed -E 's/^ +//') # Remove leading whitespace

    done <"$primer_file"
fi
echo -e "\nPrimer argument:                  $primer_arg"


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

## Run
echo -e "\n# Now running cutadapt..."
${e}Time \
    cutadapt \
        $primer_arg \
        --output "$outdir"/"$R1_basename" \
        --paired-output "$outdir"/"$R2_basename" \
        --pair-filter=any \
        $discard_arg \
        --cores "$threads" \
        $more_args \
        "$R1_in" "$R2_in"

# --pair-filter=any: Remove pair if one read is filtered (=Default)


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/"$sample_id"*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
