#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=amrfinderplus
#SBATCH --output=slurm-amrfinderplus-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "      Run AMRFinderPlus to identify antimicrobial resistance genes"
    echo "              and SNPs in a microbial genome assembly"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --nucleotide    <file>  Input nucleotide FASTA input file (genome assembly)"
    echo "  --protein       <file>  Input amino acid (protein) FASTA input file (proteome)"
    echo "  --gff           <file>  Input GFF file (annotation)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --organism      <str>   Organism name                               [default: none]"
    echo "                          For a list of options, run 'amrfinder -l'"
    echo "                          See https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option"
    echo "  --annot_format <str>    Which program produced the annotation?    [default: 'prokka']"
    echo "                          Options include 'bakta', 'patric', 'prokka', and 'rast'"
    echo "                          See https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#the---annotation_format-format-option"
    echo "  --update_db             Download/update the database"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to AMRFinderPlus"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for AMRFinderPlus and exit"
    echo "  -v/--version            Print the version of AMRFinderPlus and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --fna results/spades/smpA.fasta --aa results/prokka/smpA.faa --gff results/prokka/smpA.gff -o results/amrfinder"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/ncbi/amr/wiki"
    echo "  - Paper: https://www.nature.com/articles/s41598-021-91456-0"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/amrfinderplus
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    amrfinder --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    amrfinder --help
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
annotation_format="prokka"
update_db=false

slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly_fna=""
assembly_faa=""
gff=""
outdir=""
organism="" && organism_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )         shift && outdir=$1 ;;
        --nucleotide )          shift && assembly_fna=$1 ;;
        --protein )             shift && assembly_faa=$1 ;;
        --gff )                 shift && gff=$1 ;;
        --organism )            shift && organism=$1 ;;
        --annotation_format )   shift && annotation_format=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        --update_db )           update_db=true ;;      
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
Load_software
Set_threads

# Check input
[[ "$assembly_fna" = "" ]] && Die "Please specify an input nucleotide FASTA with --nucleotide" "$all_args"
[[ "$assembly_faa" = "" ]] && Die "Please specify an input protein FASTA with --protein" "$all_args"
[[ "$gff" = "" ]] && Die "Please specify an input GFF annotion file with --gff" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$assembly_fna" ]] && Die "Input file $assembly_fna does not exist"
[[ ! -f "$assembly_faa" ]] && Die "Input file $assembly_faa does not exist"
[[ ! -f "$gff" ]] && Die "Input file $gff does not exist"

# Define output files etc
asm_basename=$(basename "$assembly_fna")
sample_id=${asm_basename%.*}
outfile="$outdir"/"$sample_id".txt
mutation_report="$outdir"/"$sample_id"_mutation-report.txt

# Build organism argument
[[ $organism != "" ]] && organism_arg="--organism $organism"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT AMRFINDERPLUS.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input nucleotide FASTA file:      $assembly_fna"
echo "Input protein FASTA file:         $assembly_faa"
echo "Input GFF file:                   $gff"
echo "Annotation (GFF) format:          $annotation_format"
echo "Output dir:                       $outdir"
echo "Output mutation report:           $mutation_report"
echo "Update the database, if needed:   $update_db"
[[ $organism != "" ]] && echo "Organism:                         $organism"
[[ $more_args != "" ]] && echo "Other arguments for AmrFinderPlus:    $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$assembly_fna" "$assembly_faa" "$gff"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo "# Creating the output directories..."
mkdir -pv "$outdir"/logs

# Download latest database, if there's an update
if [[ "$update_db" = true ]]; then
    echo -e "\n# Checking for (& downloading) database updates...."
    Time amrfinder --update
    # Since we can't change the DB dir, this will only work if you have write access to the Conda env
fi

# Run
echo -e "\n# Running AmrfinderPlus...."
Time amrfinder \
    --nucleotide "$assembly_fna" \
    --protein "$assembly_faa" \
    --gff "$gff" \
    --annotation_format "$annotation_format" \
    --threads "$threads" \
    -o "$outfile" \
    --mutation_all "$mutation_report" \
    --name "$sample_id" \
    $organism_arg \
    $more_args


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/"$sample_id"*
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
echo
