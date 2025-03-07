#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=mpileup
#SBATCH --output=slurm-mpileup-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run 'bcftools mpileup' and 'bcftools call' to call variants using BAM files and a reference FASTA"
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA=/fs/ess/PAS0471/jelmer/conda/bcftools
readonly SCRIPT_VERSION="0.1"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY=bcftools
readonly TOOL_NAME=bcftools
readonly TOOL_DOCS=https://samtools.github.io/bcftools/bcftools.html
readonly TOOL_PAPER=https://academic.oup.com/gigascience/article/10/2/giab008/6137722
readonly VERSION_COMMAND=

# Constants - parameters
VARCALL_METHOD="--multiallelic-caller"

# Option defaults
allsites=false      # Only output variant sites

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo "  NOTE: The variant-calling method is set to '--multiallelic-caller'"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 --fasta data/assembly.fa --bam_dir results/bwa -o results/mpileup.vcf"
    echo "  - To just print the help message for this script:"
    echo "      bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--fasta      <file>  Input reference FASTA file"
    echo "  --bam_dir       <dir>   Dir with output BAM files"
    echo "  -o/--vcf_out    <file>  Output gzipped VCF file (extension 'vcf.gz')"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --allsites              Output ALL sites in the VCF, not just variants"
    echo "  --args_mpileup  <str>   Quoted string with additional argument(s) for bcftools mpileup"
    echo "  --args_call     <str>   Quoted string with additional argument(s) for bcftools call"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    local is_slurm=$1

    # Determine the location of this script, and based on that, the function script
    if [[ "$is_slurm" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/bash_functions.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
    # shellcheck source=/dev/null
    source "$function_script"
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
fasta=
bam_dir=
vcf=
args_mpileup=
args_call=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --fasta )      shift && readonly fasta=$1 ;;
        --bam_dir )         shift && readonly bam_dir=$1 ;;
        -o | --vcf )        shift && readonly vcf=$1 ;;
        --allsites )        allsites=true ;;
        --args_mpileup )    shift && readonly args_mpileup=$1 ;;
        --args_call )       shift && readonly args_call=$1 ;;
        -h )                script_help; exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$fasta" ]] && die "No input FASTA file specified, do so with -i/--fasta" "$all_args"
[[ -z "$bam_dir" ]] && die "No input BAM dir specified, do so with --bam_dir" "$all_args"
[[ -z "$vcf" ]] && die "No output VCF file specified, do so with -o/--vcf" "$all_args"
[[ ! -f "$fasta" ]] && die "Input FASTA file $fasta does not exist"
[[ ! -d "$bam_dir" ]] && die "Input BAM dir $bam_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict bash settings
set -euo pipefail

# Outputs based on script parameters
outdir=$(dirname "$vcf")
if [[ "$allsites" == true ]]; then allsites_opt=; else allsites_opt="-v"; fi

# Logging files
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output VCF file:                          $vcf"
echo "Input FASTA file:                         $fasta"
echo "Input BAM dir:                            $bam_dir"
echo "Number of BAM files:                      $(ls "$bam_dir"/*bam | wc -l)"
echo "Output all sites in the VCF:              $allsites"
[[ -n $args_mpileup ]] && echo "Other arguments for bcftools mpileup:   $args_mpileup"
[[ -n $args_call ]] && echo "Other arguments for bcftools mpileup:   $args_call"
log_time "Listing the input file(s):"
ls -lh "$fasta" "$bam_dir"/*bam
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY mpileup \
    --fasta-ref "$fasta" \
    --output-type u \
    $args_mpileup \
    "$bam_dir"/*bam |
        bcftools call \
            "$VARCALL_METHOD" \
            $allsites_opt \
            --threads "$threads" \
            --output-type z \
            $args_call \
            --output "$vcf"

#? '--output-type u/z' => uncompressed BCF / compressed VCF
#? -m: default calling methods (vs -c, old consensus calling method)

# List the output
log_time "Listing the output VCF file:"
ls -lh "$vcf"

# ==============================================================================
#                               WRAP UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
print_version "$VERSION_COMMAND" | tee "$VERSION_FILE"
script_version "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL" | tee -a "$VERSION_FILE" 
env | sort > "$ENV_FILE"
[[ "$IS_SLURM" = true ]] && resource_usage
log_time "Done with script $SCRIPT_NAME\n"
