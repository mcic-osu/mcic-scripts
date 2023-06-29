#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --out=slurm-multiqc-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "   Run MultiQC to summarize output by e.g. FastQC, Cutadapt, STAR"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <file>  Input directory (e.g. with FastQC output files)"
    echo "  -o/--outdir     <dir>   Output directory for MultiQC report"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --filename      <file>  Provide filename for MultiQC HTML report    [default: <outdir>/multiqc_report.html]"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to MultiQC"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for TODO_THIS_SOFTWARE and exit"
    echo "  -v/--version            Print the version of TODO_THIS_SOFTWARE and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/fastqc -o results/multiqc"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - '--interactive'      Reports will always have interactive figures"
    echo "  - '--force'            Overwrite existing MultiQC reports with the same name"
    echo
    echo "OUTPUT:"
    echo "  - The main output is an HTML file, by default named <outdir>/multiqc_report.html"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/multiqc
}

## Print version
Print_version() {
    Load_software
    multiqc --version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    multiqc --help
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
indir=""
outdir=""
more_args=""
filename="" && filename_arg=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --filename )        shift && filename=$1 ;;
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

## Bash script settings
set -euo pipefail

## Filename argument
[[ "$filename" != "" ]] && filename_arg="--filename $filename"

## Check input
[[ "$indir" = "" ]] && Die "Please specify an input file with -i/--indir" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"

## Report
echo
echo "=========================================================================="
echo "              STARTING SCRIPT MULTIQC.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input dir:                        $indir"
echo "Output dir:                       $outdir"
[[ $filename_arg != "" ]] && echo "Output filename:                  $filename"
[[ $more_args != "" ]] && echo "Other arguments for MultiQC:      $more_args"
echo
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

echo "## Starting MultiQC run..."
${e}Time multiqc \
    --interactive \
    --force \
    -o "$outdir" \
    $filename_arg \
    "$indir"

#? --interactive will ensure interactive plots, regardless of number of samples
#? --force will overwrite any old report


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
