#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=nfcore_ampliseq
#SBATCH --output=slurm-nfcore_ampliseq-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run the Nextflow-core metabarcoding pipeline from https://nf-co.re/ampliseq"
readonly MODULE=miniconda3
readonly CONDA=/fs/project/PAS0471/jelmer/conda/nextflow
readonly SCRIPT_VERSION="2023-07-18"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY="nextflow run"
readonly VERSION_COMMAND="nextflow -v"

# Constants - parameters
WORKFLOW_NAME=ampliseq      # The name of the nf-core workflow
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
OSC_CONFIG=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

# Parameter defaults
workflow_version=2.6.1      # The version of the nf-core workflow
workflow_dir_base=workflows/nfcore-ampliseq
workflow_dir_full=workflows/nfcore-ampliseq/${workflow_version//./_}
is_ITS=false && ITS_arg=
container_dir=/fs/project/PAS0471/containers
work_dir=/fs/scratch/PAS0471/$USER/nfcore-ampliseq
profile="singularity"
resume=true && resume_arg="-resume"

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "                $0 (v. $SCRIPT_VERSION):"
    echo "                      $DESCRIPTION"
    echo "        =============================================="
    echo "ABOUT:"
    echo "  - All workflow parameter defaults in this script are the same as in the https://nf-co.re/ampliseq pipeline,"
    echo "    except that '--ignore_failed_trimming' is always used so the pipeline will keep running when some samples have too small filesizes after trimming."
    echo "  - Different from the Nextflow default, this script will try to 'resume' (rather than restart) a previous incomplete run by default"
    echo "  - This workflow can be used for both 16S and ITS data (use the '--its' flag for the latter)"
    echo "  - Not all nf-core/ampliseq parameters are present as options to this script, use '--more_args' for parameters that aren't listed"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT \\"
    echo "    --more-args \"--extension '/*_R{1,2}.fastq.gz' --illumina_novaseq\""
    echo
    echo "NF-CORE/AMPLISEQ WORKFLOW PARAMETERS:"
    echo "Note: These are all the same as the parameters specified at https://nf-co.re/ampliseq/parameters: see there for details"
    echo "Note: To pass nf-core/ampliseq parameters that are not listed below, use the 'more_args' option."
    echo "  --input             <dir>   Input FASTQ dir                                             [REQUIRED]"
    echo "  --outdir            <dir>   Output directory (will be created)                          [REQUIRED]"
    echo "  --FW_primer         <str>   Forward primer sequence                                     [REQUIRED]"
    echo "  --RV_primer         <str>   Reverse primer sequence                                     [REQUIRED]"
    echo "  --extension         <dir>   FASTQ file extension glob pattern                           [default: '/*_R{1,2}_001.fastq.gz']"
    echo "  --trunclenf         <int>   Truncate forward reads at this length                       [default: none => uses quality instead]"
    echo "  --trunclenr         <int>   Truncate reverse reads at this length                       [default: none => uses quality instead]"
    echo "  --min_len_asv       <int>   Minimum ASV length in basepairs                             [default: none => no length filtering]"
    echo "  --max_len_asv       <int>   Maximum ASV length in basepairs                             [default: none => no length filtering]"
    echo "  --sample_inference  <str>   Dada sample inference method                                [default: 'independent']"
    echo "                                - 'independent', 'pooled', or 'pseudo'"
    echo "  --dada_ref_taxonomy <str>   Reference taxonomy for dada             [default: 'silva=138' for 16S / 'unite-fungi=8.3' for ITS]"
    echo "  --filter_ssu        <str>   Quoted, comma-separated list of kingdoms to keep after Barrnap classification"
    echo "                                - Recommended for 16S: 'bac'"
    echo "                                - Don't use this for ITS"
    echo "                                - Default: no SSU filtering"
    echo "                                - See https://nf-co.re/ampliseq/2.4.1/parameters#filter_ssu"
    echo "  --exclude_taxa      <str>   Quoted, comma-separated list of Kingdoms/taxa to remove after ASV tax. classification with DADA2 or QIIME2"
    echo "                                - Default: 'mitochondria,chloroplast'"
    echo "  --min_frequency     <int>   ASV must be present at least x times                        [default: 1]"
    echo "  --metadata          <file>  Metadata TSV file                                           [default: no metadata => no by-group stats]"
    echo "                                - At a minimum should have a column named 'id'"
    echo "                                - See https://docs.qiime2.org/2022.11/tutorials/metadata/"
    echo "  --metadata_category <str>   Metadata category/ies for downstream analysis [default: none]"
    echo "  --metadata_category_barplot <str> Metadata category/ies for barplots [default: none / 'metadata_category' if that is specified]"
    echo
    echo "OPTIONS OF THIS SHELL SCRIPT:"
    echo "  --its                       Use this flag if the data is ITS-derived, in which case:"
    echo "                                - The default taxonomy will be changed to 'unite-fungi=8.3'"
    echo "                                - The '--illumina_pe_its' and '--addsh' flags to the workflow will be used"
    echo "  --workflow_dir      <dir>   Dir with (or for) the nfcore/ampliseq workflow              [default: 'workflows/nfcore-ampliseq']"
    echo "  --workflow_version  <str>   nf-core ampliseq workflow version                           [default: $workflow_version]"
    echo "  --more_args         <str>   Additional arguments, all in a single quoted string, to pass to 'nextflow run'"
    echo 
    echo "NEXTFLOW OPTIONS:"
    echo "  -restart                    Don't attempt to resume workflow run, but start over        [default: resume workflow]"
    echo "  -config             <file>  Additional config file                                      [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  -profile            <str>   'Profile' to use from one of the config files               [default: 'singularity']"
    echo "  -work-dir           <dir>   Scratch (work) dir for the workflow                         [default: '/fs/scratch/PAS0471/\$USER/nfcore-ampliseq']"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of Nextflow and exit"
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
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
fastq_dir=
outdir=
FW_primer=
RV_primer=
exclude_taxa= && exclude_taxa_arg=
extension= && extension_arg=
trunclenf= && trunclenf_arg=
trunclenr= && trunclenr_arg=
sample_inference= && sample_inference_arg=
dada_ref_taxonomy= && dada_ref_taxonomy_arg=
min_len_asv= && min_len_asv_arg=
max_len_asv= && max_len_asv_arg=
filter_ssu= && filter_ssu_arg=
config_file=
metadata= && metadata_arg=
metadata_category=
metadata_category_barplot=
min_frequency=1
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --input )                       shift && fastq_dir=$1 ;;
        --outdir )                      shift && outdir=$1 ;;
        --extension )                   shift && extension=$1 ;;
        --FW_primer )                   shift && FW_primer=$1 ;;
        --RV_primer )                   shift && RV_primer=$1 ;;
        --its )                         is_ITS=true ;;
        --sample_inference )            shift && sample_inference=$1 ;;
        --dada_ref_taxonomy )           shift && dada_ref_taxonomy=$1 ;;
        --trunclenf )                   shift && trunclenf=$1 ;;
        --trunclenr )                   shift && trunclenr=$1 ;;
        --min_len_asv )                 shift && min_len_asv=$1 ;;
        --max_len_asv )                 shift && max_len_asv=$1 ;;
        --min_frequency )               shift && min_frequency=$1 ;;
        --filter_ssu )                  shift && filter_ssu=$1 ;;
        --exclude_taxa )                shift && exclude_taxa=$1 ;;
        --metadata )                    shift && metadata=$1 ;;
        --metadata_category )           shift && metadata_category=$1 ;;
        --metadata_category_barplot )   shift && metadata_category_barplot=$1 ;;
        --workflow_dir )                shift && workflow_dir_full=$1 ;;
        --container_dir )               shift && container_dir=$1 ;;
        --more_args )                   shift && more_args=$1 ;;
        -config )                       shift && config_file=$1 ;;
        -profile )                      shift && profile=$1 ;;
        -work-dir )                     shift && work_dir=$1 ;;
        -restart )                      resume=false && resume_arg= ;;
        -h | --help )                   script_help; exit 0 ;;
        -v )                            script_version; exit 0 ;;
        --version )                     load_env "$MODULE" "$CONDA"
                                        tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                             die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$fastq_dir" ]] && die "No FASTQ dir specified, do so with -i/--input" "$all_args"
[[ -z "$FW_primer" ]] && die "No forward primer specified, do so with --FW_primer" "$all_args"
[[ -z "$RV_primer" ]] && die "No reverse primer specified, do so with --RV_primer" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -d "$fastq_dir" ]] && die "Input dir $fastq_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
nextflow_setup
set_threads "$IS_SLURM"

# ITS options
if [[ "$is_ITS" = true ]]; then
    ITS_arg="--illumina_pe_its --addsh" #TODO declare at top
    [[ -z "$dada_ref_taxonomy" ]] && dada_ref_taxonomy='unite-fungi=8.3' #TODO declare at top
fi

# Metadata options
if [[ -n "$metadata" ]]; then
    if [[ -n "$metadata_category" ]]; then
        [[ -z "$metadata_category_barplot" ]] && metadata_category_barplot="$metadata_category"
        metadata_arg="--metadata $metadata --metadata_category $metadata_category --metadata_category_barplot $metadata_category_barplot"
    else
        metadata_arg="--metadata $metadata"
    fi
fi

# Other options - only pass these args if non-default parameters were passed
[[ -n "$sample_inference" ]] && sample_inference_arg="--sample_inference $sample_inference"
[[ -n "$dada_ref_taxonomy" ]] && dada_ref_taxonomy_arg="--dada_ref_taxonomy $dada_ref_taxonomy"
[[ -n "$min_len_asv" ]] && min_len_asv_arg="--min_len_asv $min_len_asv"
[[ -n "$max_len_asv" ]] && max_len_asv_arg="--max_len_asv $max_len_asv"
[[ -n "$trunclenf" ]] && trunclenf_arg="--trunclenf $trunclenf"
[[ -n "$trunclenr" ]] && trunclenr_arg="--trunclenr $trunclenr"
[[ -n "$min_frequency" ]] && min_frequency_arg="--min_frequency $min_frequency"
[[ -n "$filter_ssu" ]] && filter_ssu_arg="--filter_ssu $filter_ssu"
[[ -n "$extension" ]] && extension_arg="--extension $extension"
[[ -n "$exclude_taxa" ]] && exclude_taxa_arg="--exclude_taxa $exclude_taxa"

# Build the config argument
[[ ! -f "$OSC_CONFIG" ]] && OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_arg="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Other output dirs
trace_dir="$outdir"/pipeline_info

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:                 $all_args"
echo
echo "INPUT AND OUTPUT:"
echo "Input FASTQ dir:                              $fastq_dir"
echo "Output dir:                                   $outdir"
[[ -n "$extension" ]] && echo "File extension:                               $extension"
[[ -n "$metadata" ]] && echo "Metadata file:                                $metadata"
echo
echo "SETTINGS:"
echo "Forward primer:                               $FW_primer"
echo "Reverse primer:                               $RV_primer"
[[ -n "$dada_ref_taxonomy" ]] && echo "Dada reference taxonomy:                      $dada_ref_taxonomy"
[[ -n "$sample_inference" ]] && echo "Sample inference method:                      $sample_inference"
[[ -n "$min_len_asv" ]] && echo "Min. ASV length:                              $min_len_asv"
[[ -n "$max_len_asv" ]] && echo "Max. ASV length:                              $max_len_asv"
[[ -n "$trunclenf" ]] && echo "Truncate forward reads at:                    $trunclenf"
[[ -n "$trunclenr" ]] && echo "Truncate reverse reads at:                    $trunclenr"
[[ -n "$filter_ssu" ]] && echo "Filter SSU (Barrnap classification):          $filter_ssu"
[[ -n "$exclude_taxa" ]] && echo "Exclude taxa (DADA2/QIIME classification):    $exclude_taxa"
[[ -n "$min_frequency" ]] && echo "Min. ASV frequency:                           $min_frequency"
[[ -n "$metadata_category" ]] && echo "Metadata categories:                          $metadata_category"
[[ -n "$metadata_category_barplot" ]] && echo "Metadata categories for barplots:             $metadata_category_barplot"
[[ -n "$metadata_arg" ]] && echo "Metadata argument:                            $metadata_arg"
[[ -n "$more_args" ]] && echo "Additional arguments:                         $more_args"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run (if any):                 $resume"
echo "Container dir:                                $container_dir"
echo "Scratch (work) dir:                           $work_dir"
echo "Nextflow workflow dir:                        $workflow_dir_full"
echo "Config 'profile':                             $profile"
echo "Config file argument:                         $config_arg"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
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
    echo -e "\n# Downloading workflow to $workflow_dir_base"
    nf-core download "$WORKFLOW_NAME" \
        --revision "$workflow_version" \
        --compress none \
        --container-system singularity \
        --parallel-downloads "$threads" \
        --outdir "$workflow_dir_base"
    echo
fi

# Run the tool
log_time "Starting the workflow.."
runstats $TOOL_BINARY \
    "$workflow_dir_full" \
    --input "$fastq_dir" \
    $extension_arg \
    --outdir "$outdir" \
    --FW_primer "$FW_primer" \
    --RV_primer "$RV_primer" \
    --ignore_failed_trimming \
    $ITS_arg \
    $dada_ref_taxonomy_arg \
    $trunclenf_arg \
    $trunclenr_arg \
    $sample_inference_arg \
    $min_len_asv_arg \
    $max_len_asv_arg \
    $min_frequency_arg \
    $filter_ssu_arg \
    $exclude_taxa_arg \
    $metadata_arg \
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

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
