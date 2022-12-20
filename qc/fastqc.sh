#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastqc
#SBATCH --output=slurm-fastqc-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              RUN FASTQC FOR ONE OR MORE FASTQ FILES"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  Syntax: $0 -o <output-dir> [ -i <input-FASTQ> ] [...] [fastq-file-1 fastq-file-2]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  A single input FASTQ file"
    echo "                          Alternatively, pass 1 or more FASTQ files as positional arguments at the end of the command."
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to FastQC"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for FastQC and exit"
    echo "  -v/--version            Print the version of FastQC and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/A_R1.fastq -o results/fastqc"
    echo
    echo "OUTPUT:"
    echo "  - One HTML and one ZIP file per input FASTQ file"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
    echo
}


## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/fastqc-0.11.9
    set -u
}

## Print version
Print_version() {
    set +e
    Load_software
    fastqc --version
    set -e
}

## Print help for the focal program
Print_help_program() {
    Load_software
    fastqc --help
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
declare -a infiles
outdir=""
more_args=""

## Parse command-line args
all_args="$*"
count=0

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infiles[0]=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more-args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 infiles[$count]=$1 && count=$(( count + 1 )) ;;
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

## Check input
#[[ "${#infiles[@]}" -eq 0 ]] && Die "Please specify input file(s) with -i/--infile or as positional arguments" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT FASTQC.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Number of input files:            ${#infiles[@]}"
echo "Output dir:                       $outdir"
[[ $more_args != "" ]] && echo "Other arguments for FastQC:       $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
count=0
for infile in "${infiles[@]}"; do
    [[ ! -f "$infile" ]] && Die "Input file $infile does not exist"
    sampleIDs[count]=$(basename "$infile" | sed -E 's/.f?a?s?t?q.*//')
    ls -lh "$infile"
    count=$((count + 1))
done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

## Run
echo -e "\n# Now running FastQC..."
${e}Time \
    fastqc \
        --outdir "$outdir" \
        --threads "$threads" \
        "${infiles[@]}" \
        $more_args


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    for sampleID in "${sampleIDs[@]}"; do
        ls -lhd "$PWD"/"$outdir"/"$sampleID"*
    done
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
