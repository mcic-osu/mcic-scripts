#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=repeatmasker
#SBATCH --output=slurm-repeatmasker-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                        RUN REPEATMASKER"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 $0 -i <genome-FASTA> -l <genome-lib-FASTA> -o <output-dir> -s <species> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly_in   <file>   Input assembly: a nucleotide FASTA file"
    echo "  -o/--assembly_out  <file>   Output, masked, assembly"
    echo
    echo "ONE OF THESE TWO IS REQUIRED (THEY ARE MUTUALLY EXCLUSIVE):"
    echo "  --genome_lib    <file>  Genome repeat library FASTA file produced by RepeatModeler (repeatmodeler.sh script)"
    echo "  --species       <str>   Species or taxonomic group name"
    echo "                          To check which species/groups are available, run, e.g:"
    echo "                          /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'oomycetes'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to RepeatMasker"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for RepeatMasker and exit"
    echo "  -v/--version            Print the version of RepeatMasker and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/genome.fa -o results/repeatmasker"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://www.repeatmasker.org/"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    RepeatMasker --help | head -n 1
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    RepeatMasker
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
    date
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
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly_in=""
assembly_out=""
genome_lib=""
species="" && species_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly_in )    shift && assembly_in=$1 ;;
        -o | --assembly_out )   shift && assembly_out=$1 ;;
        --genome_lib )          shift && genome_lib=$1 ;;
        --species )             shift && species=$1 ;;
        --more_args )           shift && more_args=$1 ;;
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
[[ "$assembly_in" = "" ]] && Die "Please specify an input file with -i/--assembly_in" "$all_args"
[[ "$assembly_out" = "" ]] && Die "Please specify an output assembly with -o/--assembly_out" "$all_args"
[[ ! -f "$assembly_in" ]] && Die "Input file $assembly_in does not exist"
[[ "$species" = "" ]] && [[ "$genome_lib" = "" ]] && Die "Specify one of --species or --genome_lib"
[[ "$species" != "" ]] && [[ "$genome_lib" != "" ]] && Die "Specify --species or --genome_lib, not both"
[[ "$genome_lib" != "" && ! -f "$genome_lib" ]] && Die "Input file $genome_lib does not exist"

# Make file paths absolute
[[ ! "$assembly_in" =~ ^/ ]] && assembly_in=$(realpath "$assembly_in")
[[ ! "$genome_lib" =~ ^/ ]] && genome_lib=$(realpath "$genome_lib")
[[ ! "$assembly_out" =~ ^/ ]] && assembly_out="$PWD"/"$assembly_out"

# Species/genome lib args
[[ "$species" != "" ]] && species_arg="-species $species"
[[ "$genome_lib" != "" ]] && genome_lib_arg="-lib $genome_lib"

# Get outdir
outdir=$(dirname "$assembly_out")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT REPEATMASKER.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Input assembly (nucleotide FASTA):    $assembly_in"
echo "Output assembly:                      $assembly_out"
echo "Genome database from RepeatModeler:   $genome_lib"
[[ $species_arg != "" ]] && echo "Species name:                         $species"
[[ $more_args != "" ]] && echo "Other arguments for RepeatMasker:     $more_args"
echo "Number of threads/cores:              $threads"
echo
echo "Listing the input file(s):"
ls -lh "$assembly_in" "$genome_lib"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
${e}mkdir -pv "$outdir"/logs

# Move into the output dir
cd "$outdir" || exit

# Run
echo -e "\n# Now runnning RepeatMasker..."
${e}Time \
    RepeatMasker \
    -dir . \
    $genome_lib_arg \
    $species_arg \
    $more_args \
    "$assembly_in"

# Copy the output file
echo -e "\n# Now copying the masked assembly..."
cp -v "$outdir"/"$(basename "$assembly_in")".masked "$assembly_out"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee logs/version.txt
    echo -e "\n# Listing the output assembly file:"
    ls -lh "$(realpath "$assembly_out")"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo

## To check available species, e.g:
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'oomycetes'
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'stramenopiles'
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'phytophthora'
#>67593 Phytophthora megasperma f. sp. glycinea (includes), Phytophthora sojae (scientific name)
