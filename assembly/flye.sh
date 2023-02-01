#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=43
#SBATCH --mem=172G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=flye
#SBATCH --output=slurm-flye-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "           RUN FLYE TO ASSEMBLE A GENOME WITH LONG READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTQ(s)> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile    <str>   Output assembly FASTA file (use extension '.fa' or '.fasta')"
    echo "To specify the input, use one of the two following options:"
    echo "  -i/--reads      <file>  An input FASTQ file"
    echo "                          Or optionally multiple files, quoted and space-separated"
    echo "  --fofn          <file>  Text file with list of input FASTQ files one per line (fofn)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --genome-size   <str>   Genome size estimate, e.g '4.6m' or '1g'    [default: no estimate]"
    echo "  --iterations    <int>   Number of polishing iterations              [default: 1]"
    echo "  --resume                Resume previous run"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Flye"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of Flye and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/minion/my.fastq -o results/flye"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md"
    echo
}

# Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/flye-2.9.1
}

# Print version
Print_version() {
    Load_software
    flye --version
}

# Print help for the focal program
Print_help_program() {
    Load_software
    flye --help
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
    echo "Memory (MB per node): $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
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
iterations=1
resume=false && resume_arg=""

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
declare -a infiles
fofn=""
outfile=""
more_args=""
genome_size="" && genome_size_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )      shift && IFS=" " read -r -a infiles <<< "$1" ;;
        --fofn )            shift && fofn=$1 ;;
        -o | --assembly )   shift && outfile=$1 ;;
        --genome-size )     shift && genome_size=$1 ;;
        --iterations )      shift && iterations=$1 ;;
        --resume )          resume=true ;;
        --more-args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit ;;
        -h )                Print_help; exit 0;;
        --help )            Print_help_program; exit 0;;
        --debug )           debug=true ;;
        --dryrun )          dryrun=true ;;
        * )                 Print_help; Die "Invalid option $1" ;;
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

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ ${#infiles[@]} = 0 && "$fofn" = "" ]] && Die "Please specify input files with -i/--reads or --reads_fofn" "$all_args"
[[ "$outfile" = "" ]] && Die "Please specify an output file with -o/--assembly" "$all_args"

# Define output dir
outdir=$(dirname "$outfile")

# Build other args
[[ "$genome_size" != "" ]] && genome_size_arg="--genome-size $genome_size"
[[ "$resume" = true ]] && resume_arg="--resume"

# If a FOFN was provided, read file list into an array
[[ "$fofn" != "" ]] && mapfile -t infiles <"$fofn"

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FLYE.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo
[[ "$fofn" != "" ]] && echo "File with list of FASTQs (fofn):      $fofn"
echo "Input files with reads:               ${infiles[*]}"
echo "Number of input files:                ${#infiles[*]}"
echo "Output assembly file:                 $outfile"
echo "Genome size:                          $genome_size"
echo "Nr of polishing iterations:           $iterations"
echo "Resume previous run:                  $resume"
[[ $more_args != "" ]] && echo "Other arguments for Flye:             $more_args"
echo "# Listing the input files:"
for infile in "${infiles[@]}"; do
    [[ ! -f $infile ]] && Die "Input file $infile does not exist!"
    ls -lh "$infile"
done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    # Create the output directory
    echo -e "\n# Creating the output directories..."
    mkdir -pv "$outdir"/logs

    echo -e "\n# Running Flye..."
    Time flye \
        --nano-raw "${infiles[@]}" \
        --out-dir "$outdir" \
        --iterations "$iterations" \
        --threads "$threads" \
        $resume_arg \
        $genome_size_arg \
        $more_args

    # Copy assembly FASTA
    echo -e "\n# Copying the assembly FASTA file:"
    cp -v "$outdir"/assembly.fasta "$outfile"
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing the final assembly file:"
    ls -lhd "$PWD"/"$outfile"
    echo
    [[ "$slurm" = true ]] && Resource_usage
    echo
fi
echo "# Done with script"
date
