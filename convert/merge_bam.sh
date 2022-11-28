#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=merge_bam
#SBATCH --output=slurm-merge_bam-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               MERGE BAM FILES (AND SORT THE OUTPUT BAM)"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -o <output file> [ -o <input dir> / [BAM-file-1 BAM-file-2 ...]]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile    <dir>   Output BAM file"
    echo "  -i/--infile     <file>  Input dir with BAM files (all BAM files in the dir will be used)"
    echo "                          Specify the input either with an input dir or BAM files as positional arguments"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Samtools and exit"
    echo "  -v/--version            Print the version of Samtools and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  Using an input dir:"
    echo "    sbatch $0 -o results/merged.bam -i results/mapped"
    echo "  Using BAM files as positional arguments:"
    echo "    sbatch $0 -o results/merged.bam results/mapped/A.bam results/mapped/B.bam results/mapped/C.bam"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Samtools merge: http://www.htslib.org/doc/samtools-merge.html"
    echo "  - Samtools sort: http://www.htslib.org/doc/samtools-merge.html"
    echo
}

## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/samtools
    set -u
}

## Print version
Print_version() {
    Load_software
    samtools --version | head -n 2
}

## Print help for the focal program
Print_help_program() {
    Load_software
    samtools --help
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
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
declare -a bams
outdir=""
indir=""

## Parse command-line args
all_args="$*"
count=0

while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -o | --outfile )    shift && outfile=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 bams[$count]=$1 && count=$(( count + 1 )) ;;
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

## Determine file name for unsorted BAM file
outdir=$(dirname "$outfile")

## Check input
[[ "$outfile" = "" ]] && Die "Please specify an output file with -o/--outfile" "$all_args"
[[ "$indir" = "" && ${#bams[@]} = 0 ]] && Die "Please specify either an input dir with -i/--indir, or BAM files as positional args" "$all_args"
[[ "$indir" != "" && ! -d "$indir" ]] && Die "Input file $indir does not exist"

## Get BAM files if an input dir is provided
[[ "$indir" != "" ]] && mapfile bams < <(find "$indir" -type f)

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT MERGE_BAM.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Output BAM file:                  $outfile"
echo "Number of threads/cores:          $threads"
echo "Number of input BAM files:        ${#bams[@]}"
echo
echo "Listing the input BAM files:"
for bam in "${bams[@]}"; do ls -lh $bam; done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

echo "## Now running samtools..."
${e}Time samtools merge -o - ${bams[@]} |
    samtools sort -@ "$threads" > "$outfile"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/merge_bam_version.txt
    echo -e "\n# Listing the output file:"
    ls -lhd "$PWD"/"$outfile"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
