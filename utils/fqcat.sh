#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=concat_fq
#SBATCH --output=slurm-concat_fq-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                       $0"
    echo "                  CONCATENATE ONT FASTQ FILES"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input dir> -o <output file> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir     <dir>    Input file"
    echo "  -o/--outfile   <file>   Output file (its dir will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS"
    echo "  --subdir       <dir>    Subdir which may or may not be one layer removed from the indir specified by -i/--indir"
    echo "                          This can be useful if FASTQ files are inside sample-specific folders"
    echo "  --extension    <str>    Input file extension                        [default: 'fastq.gz']"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo
}

# Print SLURM job resource usage info
Resource_usage() {
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
}

# Print SLURM job requested resources
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

# Recource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Avg Mem: %t K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h"
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
# Option defaults
extension="fastq.gz"

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
indir=""
outfile=""
subdir=""

# Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -o | --outfile )        shift && outfile=$1 ;;
        --subdir )              shift && subdir=$1 ;;
        --extension )           shift && extension=$1 ;;
        -h )                    Print_help; exit 0;;
        --dryrun )              dryrun=true && e="echo ";;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash script settings
set -euo pipefail

# Check input
[[ $indir = "" ]] && Die "Please specify an input dir with -i" "$all_args"
[[ $outfile = "" ]] && Die "Please specify an output file with -o" "$all_args"
[[ ! -d $indir ]] && Die "Input dir $indir does not exist"

# Determine outdir
outdir=$(dirname "$outfile")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT CONCAT_FQ.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input dir:                        $indir"
echo "File extension:                   $extension"
[[ $subdir != "" ]] && echo "Input subdir:                     $subdir"
echo "Output file:                      $outfile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
${e}mkdir -p "$outdir"/logs

# Find the FASTQ files
if [[ $subdir != "" ]]; then
    mapfile -t fq_files < <(find "$indir"/*/"$subdir" -name "*$extension" | sort)
else
    mapfile -t fq_files < <(find "$indir" -name "*$extension" | sort)
fi

echo -e "\n# Number of FASTQ files found: ${#fq_files[@]}"
[[ ${#fq_files[@]} -eq 0 ]] && Die "No FASTQ files were found..."

# Concatenate the FASTQ files
echo -e "\n# Now concatenating the FASTQ files..."
${e}Time cat "${fq_files[@]}" > "$outfile"

# Count the number of sequences in the output file
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Now counting the number of sequences..."
    nseqs=$(zcat "$outfile" | awk '{l++}END{print l/4}')
    echo -e "# Number of sequences in output file: $nseqs"
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing the output file:"
    ls -lh "$outfile"
    echo
    [[ "$slurm" = true ]] && Resource_usage
    echo
fi
echo "# Done with script"
date
