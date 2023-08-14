#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=nfcore_rnaseq
#SBATCH --output=slurm-nfcore_rnaseq-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run the Nextflow-core RNAseq pipeline from https://nf-co.re/rnaseq
with aligner option STAR => Salmon"
SCRIPT_VERSION="2023-08-10"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY="nextflow run"
TOOL_NAME="nextflow"
VERSION_COMMAND="nextflow -v"

# Constants - parameters
WORKFLOW_NAME=rnaseq      # The name of the nf-core workflow
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
OSC_CONFIG=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

# Defaults - generic
conda_path=/fs/project/PAS0471/jelmer/conda/nextflow
container_dir=/fs/project/PAS0471/containers

# Parameter defaults
workflow_version=3.12.0                         # The version of the nf-core workflow
workflow_dir_base=workflows/nfcore-rnaseq
workflow_dir_full="$workflow_dir_base"/${workflow_version//./_}
work_dir=/fs/scratch/PAS0471/$USER/nfcore-rnaseq
profile="singularity"
resume=true && resume_arg="-resume"

# ==============================================================================
#                                FUNCTIONS
# ==============================================================================
Print_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLES:"
    echo "  sbatch $0 -i samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/nfc_rnaseq"
    echo "  sbatch $0 -i samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/nfc_rnaseq \\"
    echo "    -a \"--skip_bbsplit false --bbsplit_fasta_list contaminant_refs.tsv\""
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--samplesheet    <file>  Sample sheet containing paths to FASTQ files and sample info"
    echo "                                - See example code below, and https://nf-co.re/rnaseq/docs/usage#samplesheet-input for more details"
    echo "  -o/--outdir         <dir>   Output directory (will be created if needed)"
    echo "  --ref_fasta         <file>  Reference genome FASTA file"
    echo "  --ref_annot         <file>  Reference genome annotation file (GTF/GFF - GTF format preferred)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --workflow_version  <str>   Nf-core rnaseq workflow version         [default: $workflow_version]"
    echo "  --workflow_dir      <dir>   Dir with/for the workflow repo          [default: $workflow_dir_base]"
    echo "                                - If the correct version of the workflow is already present in this dir, it won't be downloaded again"
    echo "  --opts              <str>   Additional options to pass to $TOOL_NAME run"
    echo
    echo "NEXTFLOW-RELATED OPTIONS:"
    echo "  --restart                   Restart workflow from the beginning     [default: resume workflow if possible]"
    echo "  --container_dir     <dir>   Directory with container images         [default: $container_dir]"
    echo "                                - Required images will be downloaded here when not already present here" 
    echo "  --config            <file>  Additional config file                  [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --profile            <str>  'Profile' to use from one of the config files [default: $profile]"
    echo "  --work_dir           <dir>  Scratch (work) dir for the workflow     [default: $work_dir]"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of Nextflow and exit"
    echo
    echo "HARDCODED WORKFLOW OPTIONS:"
    echo "  --aligner star_salmon"
    echo "  --remove_ribo_rna"
    echo "  --save_reference"
    echo "  --save_non_ribo_reads"
    echo "  --save_merged_fastq"
    echo
    echo "SOME KEY OUTPUT FILES:"
    echo "  - HTML file with summary of results: <outdir>/multiqc/star_salmon/multiqc_report.html"
    echo "  - Gene counts for use with DESeq2 in <outdir>/star_salmon/salmon.merged.gene_counts_length_scaled.rds"
    echo "    Use the script 'mcic-scripts/R_templates/nfcore-rnaseq_load-counts.R' to get started with that."
    echo
    echo "NFCORE RNASEQ WORKFLOW DOCUMENTATION:"
    echo "   - https://nf-co.re/rnaseq "
    echo
    echo "EXAMPLE CODE FOR CREATING THE SAMPLE SHEET"
    echo "  wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py"
    echo "  fqdir=data/fastq                         # Dir with your FASTQ files"
    echo "  samplesheet=data/meta/samplesheet.csv    # Output samplehseet"
    echo "  python3 fastq_dir_to_samplesheet.py \$fqdir \$samplesheet --strandedness reverse"
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
    function_script=$(realpath "$script_dir"/../dev/bash_functions2.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions2.sh
    fi
    source "$function_script"
}

nextflow_setup() {
    # Singularity container dir - any downloaded containers will be stored here;
    # if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"

    # Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
samplesheet=""
ref_annot=""
ref_fasta=""
outdir=""
config_file=""
opts=""
version_only=false
threads=

# Parse command-line options
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --samplesheet )        shift && samplesheet=$1 ;;
        -o | --outdir )             shift && outdir=$1 ;;
        --ref_annot )               shift && ref_annot=$1 ;;
        --ref_fasta )               shift && ref_fasta=$1 ;;
        --container_dir )           shift && container_dir=$1 ;;
        --opts )                    shift && opts=$1 ;;
        --config | -config )        shift && config_file=$1 ;;
        --profile | -profile )      shift && profile=$1 ;;
        --work_dir | -work-dir )    shift && work_dir=$1 ;;
        --restart | -restart )      resume=false && resume_arg= ;;
        -h | --help )               script_help; exit 0 ;;
        -v )                        script_version; exit 0 ;;
        --version )                 version_only=true ;;
        * )                         script_help; die "Invalid option $1";;
    esac
    shift
done

# Check input
[[ "$samplesheet" = "" ]] && die "Please specify a samplesheet with -i" "$all_args"
[[ "$ref_fasta" = "" ]] && die "Please specify a genome FASTA file with -f" "$all_args"
[[ "$ref_annot" = "" ]] && die "Please specify an annotation (GFF/GTF) file with -g" "$all_args"
[[ "$outdir" = "" ]] && die "Please specify an output dir with -o" "$all_args"
[[ ! -f "$samplesheet" ]] && die "Samplesheet $samplesheet does not exist"
[[ ! -f "$ref_fasta" ]] && die "Reference FASTA file $ref_fasta does not exist"
[[ ! -f "$ref_annot" ]] && die "Reference annotation file $ref_annot does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Bash strict settings
set -euo pipefail

# Load software
load_env "$conda_path"
nextflow_setup
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Build the config argument
[[ ! -f "$OSC_CONFIG" ]] && OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_arg="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Setup Nextflow arguments: annotation filetype
if [[ "$ref_annot" =~ .*\.gff3? ]]; then
    annot_arg="--gff $ref_annot"
elif [[ "$ref_annot" =~ .*\.gtf ]]; then
    annot_arg="--gtf $ref_annot"
else
    die "Unknown annotation file format"
fi

# Other output dirs
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
trace_dir="$outdir"/pipeline_info

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script: $all_opts"
echo
echo "INPUT AND OUTPUT:"
echo "Sample sheet:                     $samplesheet"
echo "Reference genome FASTA file:      $ref_fasta"
echo "Reference genome annotation file: $ref_annot"
echo "Output dir:                       $outdir"
echo
echo "SETTINGS:"
[[ -n "$opts" ]] && echo "Additional options:               $opts"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run:              $resume"
echo "Container dir:                    $container_dir"
echo "Scratch (work) dir:               $work_dir"
echo "Nextflow workflow dir:            $workflow_dir_full"
echo "Config 'profile':                 $profile"
echo "Config file argument:             $config_arg"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
set_threads "$IS_SLURM"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                              RUN
# ==============================================================================
# Make necessary dirs
log_time "Creating the output directories..."
mkdir -pv "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

# Download the OSC config file
if [[ ! -f "$OSC_CONFIG" ]]; then
    log_time "Downloading the mcic-scripts Nextflow OSC config file to $OSC_CONFIG..."
    wget -q -O "$OSC_CONFIG" "$OSC_CONFIG_URL"
fi

# Download workflow, if needed
if [[ ! -d "$workflow_dir_full" ]]; then
    mkdir -p "$(dirname "$workflow_dir_base")"
    log_time "Downloading the workflow repository to $workflow_dir_base"
    nf-core download "$WORKFLOW_NAME" \
        --revision "$workflow_version" \
        --compress none \
        --container-system singularity \
        --parallel-downloads "$threads" \
        --outdir "$workflow_dir_base"
    echo
fi

# Run the workflow
log_time "Starting the workflow.."
runstats $TOOL_BINARY \
    "$workflow_dir_full" \
    --input "$samplesheet" \
    --outdir "$outdir" \
    --fasta "$ref_fasta" \
    $annot_arg \
    --aligner star_salmon \
    --remove_ribo_rna \
    --save_reference \
    --save_non_ribo_reads \
    --save_merged_fastq \
    -work-dir "$work_dir" \
    -with-report "$trace_dir"/report.html \
    -with-trace "$trace_dir"/trace.txt \
    -with-timeline "$trace_dir"/timeline.html \
    -with-dag "$trace_dir"/dag.png \
    -ansi-log false \
    -profile "$profile" \
    $config_arg \
    $resume_arg \
    $opts

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
