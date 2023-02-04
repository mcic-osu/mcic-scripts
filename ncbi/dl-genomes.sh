#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=dl_genomes
#SBATCH --output=slurm-dl_genomes-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               DOWNLOAD GENOMES WITH THE NCBI DATASETS TOOL"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--accessions     <file>  Input file with Assembly or BioProject Accessions"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --include           <str>   'genome', 'rna', 'protein', 'cds', 'gff3', 'gtf', 'gbff', and/or 'seq-report' (comma-separated) [default: 'genome']"
    echo "  --assembly_version  <str>   'latest' or 'all'"
    echo "  --more_args         <str>   Quoted string with additional argument(s) to pass to 'datasets download genome'"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for 'datasets download genome' and exit"
    echo "  -v/--version                Print the version of datasets and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i accessions.txt -o results/ncbi_genomes"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/ncbi-datasets

    export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    datasets --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    datasets download genome --help
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
include=genome
assembly_version=latest

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
accessions=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --accessions )     shift && accessions=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        --include )             shift && include=$1 ;;
        --assembly_version )    shift && assembly_version=$1 ;;
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

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software

# Make input path absolute
accessions=$(realpath "$accessions")

# Check input
[[ "$accessions" = "" ]] && Die "Please specify an input file with -i/--accessions" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$accessions" ]] && Die "Input file $accessions does not exist"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT DL-GENOMES.SH"
date
echo "=========================================================================="
echo "All arguments to this script:                     $all_args"
echo "Output dir:                                       $outdir"
echo "Input file with accession list:                   $accessions"
echo "Number of input accessions:                       $(wc -l < "$accessions")"
[[ $more_args != "" ]] && echo "Other arguments for 'datasets download genome':   $more_args"
echo
echo "Listing the input file(s):"
ls -lh "$accessions"
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
cd "$outdir" || exit 1

# Run
echo -e "\n# Running NCBI datasets..."
${e}Time \
    datasets download genome \
        accession \
        --inputfile "$accessions" \
        --include "$include" \
        --assembly-version "$assembly_version" \
        --filename genomes.zip \
        --api-key "$NCBI_API_KEY" \
        $more_args

# Process the output
echo -e "\n# Unzipping the downloaded archive..."
unzip -o genomes.zip

echo -e "\n# Moving all files into the outdir..."
find . -type f -exec mv -v {} . \;

# Report
echo
echo "Number of input accessions:                       $(wc -l < "$accessions")"
echo "Number of output genomes ('*fna' files):          $(find . -name "*fna" | wc -l)"


#? Alternative:
# Use AstroBioMike's bit - https://github.com/AstrobioMike/bit
#micromamba activate /fs/ess/PAS0471/jelmer/conda/bit
#bit-dl-ncbi-assemblies -w "$acc_file" -f fasta -j 1


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
