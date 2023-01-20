#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=medaka
#SBATCH --output=slurm-medaka-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              RUN MEDAKA TO POLISH A GENOME ASSEMBLY"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --reads <FASTQ> --assembly <FASTA> --model <Medaka-model> -o <output-file> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile    <file>  Output assembly FASTA (dir will be created if needed)"
    echo "  --reads         <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  --assembly      <file>  Input assembly: FASTA file (to be corrected)"
    echo "  --model         <str>   Medaka model, see the Medaka docs at https://github.com/nanoporetech/medaka#models"
    echo "                          Get a full list of possible models by running:"
    echo "                              module load miniconda3/4.12.0-py39"
    echo "                              source activate /fs/ess/PAS0471/jelmer/conda/medaka-1.7.2"
    echo "                              medaka tools list_models"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Medaka"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of Medaka and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/my.fastq -r results/assembly.fasta -o results/medaka -m r941_min_hac_g507"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - https://github.com/nanoporetech/medaka"
    echo
}

# Load software
Load_software() {
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate 2>/dev/null; done
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/medaka-1.7.2
}

# Print version
Print_version() {
    Load_software
    medaka_consensus -h &> medaka_help.txt
    sed -n '2p' medaka_help.txt
    rm medaka_help.txt
}

# Print help for the focal program
Print_help_program() {
    Load_software
    medaka_consensus -h
}


# Print SLURM job resource usage info
Resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
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
debug=false
dryrun=false
slurm=true

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
outfile=""
reads=""
assembly=""
model=""
more_args=""

# Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -o | --outfile )   shift && outfile=$1 ;;
        --reads )          shift && reads=$1 ;;
        --assembly )       shift && assembly=$1 ;;
        --model )          shift && model=$1 ;;
        --more_args )      shift && more_args=$1 ;;
        -h )               Print_help; exit 0;;
        --help )           Print_help_program; exit 0;;
        --debug )          debug=true ;;
        --dryrun )         dryrun=true ;;
        * )                Die "Invalid option $1" "$all_args" ;;
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

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Bash script settings
set -euo pipefail

# Check input
[[ $reads = "" ]] && Die "Please specify a file with input reads with --reads" "$all_args"
[[ $assembly = "" ]] && Die "Please specify an input assembly with --assembly" "$all_args"
[[ $outfile = "" ]] && Die "Please specify an output file with -o" "$all_args"
[[ $model = "" ]] && Die "Please specify a Medaka model with --model" "$all_args"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"
[[ ! -f $assembly ]] && Die "Input assembly file $assembly does not exist"

# Get output dir
outdir=$(dirname "$outfile")

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT MEDAKA.SH"
date
echo "=========================================================================="
echo "Input reads (FASTQ) file:         $reads"
echo "Input assembly (FASTA) file:      $assembly"
echo "Medaka model:                     $model"
echo "Output assembly file:             $outfile"
[[ $more_args != "" ]] && echo "Other arguments for Medaka:    $more_args"
echo "# Listing the input files:"
ls -lh "$reads" "$assembly" 
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

    # Run
    echo -e "\n# Running Medaka..."
    Time medaka_consensus \
        -i "$reads" \
        -d "$assembly" \
        -o "$outdir" \
        -t "$threads" \
        -m "$model"

    echo -e "\n# Renaming the output file"
    mv -v "$outdir"/consensus.fasta "$outfile"
fi

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
    echo
    [[ "$slurm" = true ]] && Resource_usage
    echo
fi
echo "# Done with script"
date
