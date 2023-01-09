#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=nfcore_rnaseq
#SBATCH --output=slurm-nfcore_rnaseq-%j.out

# ==============================================================================
#                                FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================================="
    echo "                             $0"
    echo "        Run the Nextflow-core RNAseq pipeline from https://nf-co.re/rnaseq"
    echo "======================================================================================="
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <samplesheet> -f <reference-fasta> -g <reference-annot> -o <outdir> [...]"
    echo "    (Always submit the script to the queue with 'sbatch' when wanting to start a workflow run)"
    echo "  bash $0 -h"
    echo "    (Run the script using 'bash' and the -h option to get help)"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--samplesheet <file>    Sample sheet containing paths to FASTQ files and sample info"
    echo "                             (See below and https://nf-co.re/rnaseq/3.9/usage#samplesheet-input)"
    echo "  -o/--outdir     <dir>      Output directory (will be created if needed)"
    echo "  --ref_fasta     <file>     Reference genome FASTA file"
    echo "  --ref_annot     <file>     Reference genome annotation file (GTF/GFF - GTF format preferred)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args     <str>       Additional arguments to pass to 'nextflow run'"
    echo
    echo "NEXTFLOW-RELATED OPTIONS:"
    echo "  --nf_file       <file>      Main .nf workflow definition file                         [default: 'workflows/nfcore-rnaseq/workflow/main.nf']"
    echo "  -no-resume                  Don't attempt to resume workflow run, but start over      [default: resume workflow]"
    echo "  -profile        <str>       Profile to use from one of the config files               [default: 'singularity']"
    echo "  -work-dir       <dir>       Scratch (work) dir for the workflow                       [default: '/fs/scratch/PAS0471/\$USER/nfc_rnaseq']"
    echo "                                - This is where the workflow results will be stored"
    echo "                                  before final results are copied to the output dir."
    echo "  -c/-config      <file      Additional config file                                    [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --container_dir <dir>       Singularity container dir                                 [default: '/fs/project/PAS0471/containers']"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                     Turn on debugging mode: print all commands"
    echo "  -h/--help                   Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/nfc_rnaseq"
    echo "  sbatch $0 -i samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/nfc_rnaseq \\"
    echo "    -a \"--skip_bbsplit false --bbsplit_fasta_list contaminant_refs.tsv\""
    echo
    echo "HARDCODED WORKFLOW OPTIONS:"
    echo "  --aligner star_salmon"
    echo "  --remove_ribo_rna"
    echo "  --save_reference"
    echo "  --save_non_ribo_reads"
    echo "  --save_merged_fastq"
    echo "  -ansi-log false"
    echo
    echo "OUTPUT:"
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
    echo '  python3 fastq_dir_to_samplesheet.py "$fqdir" "$samplesheet" --strandedness reverse'
}

Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/nextflow
    set -u

    # Singularity container dir - any downloaded containers will be stored here;
    # if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"

    # Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
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
    echo -e "\nFor help, run this script with the '-h' / '--help' option"
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
#                     CONSTANTS AND DEFAULTS
# ==============================================================================
# URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

# Option defaults
nf_file=workflows/nfcore-rnaseq/workflow/main.nf
container_dir=/fs/project/PAS0471/containers
work_dir=/fs/scratch/PAS0471/$USER/nfc_rnaseq
profile="singularity"
resume=true && resume_arg="-resume"
debug=false


# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
samplesheet=""
ref_annot=""
ref_fasta=""
outdir=""
config_file=""
more_args=""

# Parse command-line options
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --samplesheet )    shift && samplesheet=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --ref_annot )           shift && ref_annot=$1 ;;
        --ref_fasta )           shift && ref_fasta=$1 ;;
        --nf_file )             shift && nf_file=$1 ;;
        --container_dir )       shift && container_dir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        -config )               shift && config_file=$1 ;;
        -profile )              shift && profile=$1 ;;
        -work-dir )             shift && work_dir=$1 ;;
        -no-resume )            resume=false ;;
        --debug )               debug=true ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_help; Die "Invalid option $1";;
    esac
    shift
done


# ==============================================================================
#                              OTHER SETUP
# ==============================================================================
# Bash strict settings
set -ueo pipefail

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Load Conda environment
Load_software

# Check input
[[ "$samplesheet" = "" ]] && Die "Please specify a samplesheet with -i" "$all_args"
[[ "$ref_fasta" = "" ]] && Die "Please specify a genome FASTA file with -f" "$all_args"
[[ "$ref_annot" = "" ]] && Die "Please specify an annotation (GFF/GTF) file with -g" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o" "$all_args"
[[ ! -f "$samplesheet" ]] && Die "Samplesheet $samplesheet does not exist"
[[ ! -f "$ref_fasta" ]] && Die "Reference FASTA file $ref_fasta does not exist"
[[ ! -f "$ref_annot" ]] && Die "Reference annotation file $ref_annot does not exist"

# Build the config argument
[[ ! -f "$osc_config" ]] && osc_config="$outdir"/$(basename "$OSC_CONFIG_URL")
config_arg="-c $osc_config"
[[ "$config_file" != "" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Setup Nextflow arguments: resume option
[[ "$resume" = false ]] && resume_arg=""

# Setup Nextflow arguments: annotation filetype
if [[ "$ref_annot" =~ .*\.gff3? ]]; then
    annot_arg="--gff $ref_annot"
elif [[ "$ref_annot" =~ .*\.gtf ]]; then
    annot_arg="--gtf $ref_annot"
else
    Die "Unknown annotation file format"
fi

# Other output dirs
trace_dir="$outdir"/pipeline_info

# Report
echo
echo "=========================================================================="
echo "              STARTING SCRIPT NCFORE_RNASEQ.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Sample sheet:                     $samplesheet"
echo "Reference genome FASTA file:      $ref_fasta"
echo "Reference genome annotation file: $ref_annot"
echo "Output dir:                       $outdir"
echo "Resume previous run:              $resume"
echo
echo "Container dir:                    $container_dir"
echo "Scratch (work) dir:               $work_dir"
echo "Nextflow workflow definition file:$nf_file"
echo "Config 'profile':                 $profile"
echo "Config file argument:             $config_arg"
[[ "$config_file" != "" ]] && echo "Additional config file:           $config_file"
[[ "$more_args" != "" ]] && echo "Additional arguments:             $more_args"
echo "=========================================================================="
echo

# ==============================================================================
#                              RUN
# ==============================================================================
# Make necessary dirs
echo "# Creating the output directories..."
mkdir -pv "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

# Remove old trace files
[[ -f "$trace_dir"/report.html ]] && rm "$trace_dir"/report.html
[[ -f "$trace_dir"/trace.txt ]] && rm "$trace_dir"/trace.txt
[[ -f "$trace_dir"/timeline.html ]] && rm "$trace_dir"/timeline.html
[[ -f "$trace_dir"/dag.png ]] && rm "$trace_dir"/dag.png

# Download the OSC config file
[[ ! -f "$osc_config" ]] && wget -q -O "$osc_config" "$OSC_CONFIG_URL"

# Download workflow, if needed
if [[ ! -f "$nf_file" ]]; then
    workflow_dir="$(dirname "$(dirname "$nf_file")")"
    mkdir -p "$(dirname "$workflow_dir")"
    echo "# Downloading workflow to $workflow_dir"
    nf-core download rnaseq \
        --revision 3.9 \
        --compress none \
        --container singularity \
        --outdir "$workflow_dir"
    echo
fi

# Define the workflow command
echo -e "# Starting the workflow...\n"
Time nextflow run \
    "$nf_file" \
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
    $more_args


# ==============================================================================
#                              WRAP UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "\n# Listing files in the output dir:"
ls -lhd "$PWD"/"$outdir"/*
Resource_usage
echo "# Done with script"
date
