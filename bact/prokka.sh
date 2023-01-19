#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=prokka
#SBATCH --output=slurm-prokka-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         ANNOTATE A PROKARYOTIC GENOME ASSEMBLY WITH PROKKA"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <assembly-FASTA> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly   <file>  Input genome assembly FASTA file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --genus         <str>   Genus name of the focal organism"
    echo "  --species       <str>   Quoted string with additional argument(s) to pass to Prokka"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Prokka"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Prokka and exit"
    echo "  -v/--version            Print the version of Prokka and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/spades/assembly.fasta -o results/prokka"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/tseemann/prokka"
    echo "  - Paper: https://pubmed.ncbi.nlm.nih.gov/24642063/"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/prokka-1.14.6
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    #TODO_THIS_SOFTWARE --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    #TODO_THIS_SOFTWARE --help
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
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly=""
outdir=""
genus="" && species_arg=""
species="" && genus_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly )   shift && assembly=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --genus )           shift && genus=$1 ;;
        --species )         shift && species=$1 ;;
        --more_args )       shift && more_args=$1 ;;
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
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Get the sample ID
file_ext=$(basename "$assembly" | sed -E 's/.*(.fasta|.fa|.fna)$/\1/')
sampleID=$(basename "$assembly" "$file_ext")

# Build the 'genus' and 'species' name arguments
[[ "$genus" != "" ]] && genus_arg="--genus $genus"
[[ "$species" != "" ]] && species_arg="--species $species"

# Check input
[[ "$assembly" = "" ]] && Die "Please specify an input file with -i/--assembly" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$assembly" ]] && Die "Input file $assembly does not exist"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT PROKKA.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input assembly file:              $assembly"
echo "Output dir:                       $outdir"
echo "Genus:                            $genus"
echo "Species:                          $species"
echo "Sample ID:                        $sampleID"
[[ $more_args != "" ]] && echo "Other arguments for Prokka:       $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$assembly"
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

# Run
echo -e "\n# Now running Prokka..."
${e}Time prokka \
    --outdir "$outdir" \
    --prefix "$sampleID" \
    --cpus "$threads"  \
    --force \
    $genus_arg \
    $species_arg \
    $more_args \
    "$infile"

#? --strain
#? --usegenus  ?
#? --addgenes

## Remove DNA sequences from GFF file
echo -e "\n# Now creating a copy of the GFF file without sequences..."
${e}mv -v "$outdir"/"$sampleID".gff "$outdir"/"$sampleID"_withseqs.gff
${e}sed '/^##FASTA/Q' "$outdir"/"$sampleID"_withseqs.gff > "$outdir"/"$sampleID".gff

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/"$sampleID"*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
