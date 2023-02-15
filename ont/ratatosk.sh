#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=60G
#SBATCH --job-name=ratatosk
#SBATCH --output=slurm-ratatosk-%j.out

#? This script processes about 28k ONT reads per hour

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         RUN RATATOSK TO CORRECT LONG READS WITH ILLUMINA READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <longread-FASTQ> -I <shortread-FASTQ-list> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--fq_long       <file>   Input long-read FASTQ file"
    echo "  -I/--fq_short_list <file>   List with input short-read FASTQ file(s)"
    echo "  -o/--outdir        <dir>    Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --insert_size      <int>    Insert size - set to read length for single-end reads  [default: 500]"
    echo "  --more_args        <str>    Quoted string with additional argument(s) to pass to Ratatosk"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v/--version                Print the version of Ratatosk and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/minion/my.fastq.gz -I data/illumina/fqlist.txt -o results/ratatosk"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Repo/docs: https://github.com/DecodeGenetics/Ratatosk"
    echo "  - Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02244-4"
    echo
}

# Load software
Load_software() {
    RATATOSK=/fs/ess/PAS0471/jelmer/software/ratatosk/Ratatosk
}

# Print version
Print_version() {
    Load_software
    "$RATATOSK" --version
}

# Set the number of threads/CPUs
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

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
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

# Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
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
# Option defaults
insert_size=500

slurm=true
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
fq_long=""
fq_short_list=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --fq_long )        shift && fq_long=$1 ;;
        -I | --fq_short_list )  shift && fq_short_list=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --insert_size )         shift && insert_size=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true ;;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_help; Die "Invalid option $1" ;;
    esac
    shift
done

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Load software
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Bash script settings
set -euo pipefail

# Check input
[[ $fq_long = "" ]] && Die "Please specify a long-read FASTQ file with --fq_long"
[[ $fq_short_list = "" ]] && Die "Please specify a list with short-read FASTQ files with --fq_short_list"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $fq_long ]] && Die "Input file $fq_long does not exist"
[[ ! -f $fq_short_list ]] && Die "Input file $fq_short_list does not exist"

# Define output files (NOTE: Ratatosk will add .fastq)
fq_out="$outdir"/$(basename "$fq_long" .fastq.gz)

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RATATOSK.SH"
date
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input long-read FASTQ file:               $fq_long"
echo "List with input short-read FASTQ file(s): $fq_short_list"
echo "Output dir:                               $outdir"
echo "Insert size:                              $insert_size"
echo "Output FASTQ file:                        $fq_out"
[[ $more_args != "" ]] && echo "Other arguments for Ratatosk:    $more_args"
echo "# Listing the input files:"
ls -lh "$fq_long"
echo
mapfile -t infiles_short <"$fq_short_list"
ls -lh "${infiles_short[@]}"
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    
    # Create the output directory
    mkdir -pv "$outdir"/logs

    # Run
    echo -e "\n# Now running Ratatosk..."
    Time "$RATATOSK" correct \
        -v \
        --cores "$threads" \
        --insert-sz "$insert_size" \
        --in-short "$fq_short_list" \
        --in-long "$fq_long" \
        --out-long "$fq_out"

    # Gzip FASTQ file
    echo -e "\n# Now gzipping the output FASTQ file..."
    gzip "$fq_out".fastq
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt

    echo -e "\n# Number of reads in the input file:"
    zcat "$fq_long" | awk '{ s++ } END{ print s/4 }'
    echo -e "\n# Number of reads in the output file:"
    zcat "$fq_out".fastq.gz | awk '{ s++ } END{ print s/4 }'

    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    echo
    [[ "$slurm" = true ]] && Resource_usage
fi
echo
echo "# Done with script"
date
