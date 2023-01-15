#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=longstitch
#SBATCH --output=slurm-longstitch-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "    RUN LONGSTITCH TO SCAFFOLD A GENOME ASSEMBLY WITH LONG READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --assembly <assembly FASTA> --fastq <FASTQ file> --genome_size <genome size> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --assembly      <file>  Input assembly FASTA file"
    echo "  --fastq         <file>  Input FASTQ file with long reads"
    echo "  --genome_size   <int>   Estimated genome size"
    echo "  -o/--outfile    <file>  Output assembly FASTQ file (dir will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --run_mode      <str>   Run mode: 'tigmint-ntLink-arks', 'tigmint-ntLink', or 'ntLink-arks' [default: 'tigmint-ntLink-arks']"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Longstitch"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i TODO -o results/TODO"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/bcgsc/LongStitch"
    echo "  - Paper: "
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/longstitch-1.0.3
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
    echo "Memory (MB per node): $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):      $SLURM_CPUS_PER_TASK"
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
run_mode="tigmint-ntLink-arks"

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly=""
fastq=""
genome_size=""
outfile=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --assembly )        shift && assembly=$1 ;;
        --fastq )           shift && fastq=$1 ;;
        --genome_size )     shift && genome_size=$1 ;;
        --run_mode )        shift && run_mode=$1 ;;
        -o | --outfile )    shift && outfile=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -h )                Print_help; exit 0 ;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
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

# Make paths absolute, because we have to move into the output dir
[[ ! "$assembly" =~ ^/ ]] && assembly="$PWD"/"$assembly"
[[ ! "$fastq" =~ ^/ ]] && fastq="$PWD"/"$fastq"

# Check input
[[ "$assembly" = "" ]] && Die "Please specify an input assembly FASTA with --assembly" "$all_args"
[[ "$fastq" = "" ]] && Die "Please specify an input FASTQ file with --fastq" "$all_args"
[[ "$genome_size" = "" ]] && Die "Please specify a genome size with --genome_size" "$all_args"
[[ "$outfile" = "" ]] && Die "Please specify an output file with -o/--outfile" "$all_args"
[[ ! -f "$assembly" ]] && Die "Input assembly file $assembly does not exist"
[[ ! -f "$fastq" ]] && Die "Input FASTQ file $fastq does not exist"

# Determine output prefix
outdir=$(dirname "$outfile")
file_ext=$(basename "$assembly" | sed -E 's/.*(.fasta|.fa|.fna)$/\1/')
out_prefix=$(basename "$assembly" "$file_ext")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT LONGSTITCH.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input assembly FASTA file:        $assembly"
echo "Input FASTQ file:                 $fastq"
echo "Output assembly:                  $outfile"
echo "Genome size:                      $genome_size"
[[ $more_args != "" ]] && echo "Other arguments for Longstitch:   $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$assembly" "$fastq"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo "Now creating the output directories..."
${e}mkdir -pv "$outdir"/logs

# Make links to input files
assembly_base=assembly
fastq_base=reads
ln -s "$assembly" "$outdir"/"$assembly_base".fa
ln -s "$fastq" "$outdir"/"$fastq_base".fq.gz

# Run
echo -e "\n# Now running Longstitch..."
cd "$outdir" || exit 1
${e}Time \
    longstitch \
    $run_mode \
    draft="$assembly_base" \
    reads="$fastq_base" \
    G=$genome_size \
    out_prefix="$out_prefix" \
    t="$threads" \
    $more_args

# Copy the output file
echo -e "\n# Now copying the output file..."
cp -v "$(readlink -f "$out_prefix".scaffolds.fa)" "$outfile"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
