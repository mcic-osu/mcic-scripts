#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=rcorrfilter
#SBATCH --output=slurm-rcorrfilter-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "        Filter an rcorrector-processed pair of FASTQ files"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input R1 FASTQ> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  R1 input FASTQ file as output by rcorrector.sh (name of R2 file will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/rcorrector/A_R1.fastq.gz -o results/rcorrfilter"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    #CONDA_ENV=/users/PAS0471/jelmer/miniconda3/envs/rcorrector-env
    CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/rcorrector-1.0.5
    FILTER_SCRIPT=$CONDA_ENV/bin/FilterUncorrectabledPEfastq.py
    source activate "$CONDA_ENV"
    
    #? The script `FilterUncorrectabledPEfastq.py` from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
    #? has been added to the rcorrector Conda environment
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
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
debug=false
dryrun=false && e=""
slurm=true

## Placeholder defaults
R1=""
outdir=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        -h | --help )       Print_help; exit 0 ;;
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

## Bash script settings
set -euo pipefail

## FASTQ filename parsing TODO_edit_or_remove
file_ext=$(basename "$R1" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2=${R1/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")

## Define output files
R1_out="$outdir"/$(basename "$R1" .cor.fq.gz).fastq
R2_out="$outdir"/$(basename "$R2" .cor.fq.gz).fastq
outdir_logs="$outdir"/logs

## Check input
[[ "$R1" = "" ]] && Die "Please specify an R1 input file with -i/--R1" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$R1" ]] && Die "R1 Input file $R1 does not exist"
[[ ! -f "$R2" ]] && Die "R2 Input file $R2 does not exist"
[[ "$R1" = "$R2" ]] && Die "Input file R1 and R2 refer to the same path ($R1)"

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TODO_SCRIPTNAME"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "R1 input file:                    $R1"
echo "Output dir:                       $outdir"
echo
echo "R2 input file:                    $R2"
echo "R1 output file:                   $R1_out"
echo "R2 output file:                   $R2_out"
echo
echo "Listing the input file(s):"
ls -lh "$R1" "$R2"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir" "$outdir_logs"

    ## Run the Python script
    echo "## Running read-filter script..."
    Time python2 "$FILTER_SCRIPT" -s "$sample_id" -1 "$R1" -2 "$R2"

    ## Move output files
    echo -e "\n## Moving and gzipping output files..."
    echo unfixrm*"$sample_id"*R1*cor.fq
    mv -v unfixrm_"$sample_id"*R1*cor.fq "$R1_out"
    Time gzip -fv "$R1_out"

    echo unfixrm*"$sample_id"*R2*cor.fq
    mv -v unfixrm_"$sample_id"*R2*cor.fq "$R2_out"
    Time gzip -fv "$R2_out"

    mv -v rmunfixable_"$sample_id".log "$outdir_logs"
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
