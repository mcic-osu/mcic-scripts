#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=nfc_ampliseq
#SBATCH --output=slurm-nfc_ampliseq-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run the Nextflow-core metabarcoding pipeline from https://nf-co.re/ampliseq"
SCRIPT_VERSION="2024-02-16"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY="nextflow run"
export TOOL_NAME="nextflow"
VERSION_COMMAND="nextflow -v"

# Constants - parameters
WORKFLOW_NAME=ampliseq                                  # The name of the nf-core workflow
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config

# Parameter defaults - workflow
workflow_version=2.8.0                                  # The version of the nf-core workflow
ITS_taxonomy='unite-fungi=8.3'
is_ITS=false && ITS_opt=
ITS_opt_default="--illumina_pe_its --addsh"             # When is_ITS is true, use this arg

# Parameter defaults - infrastructure
osc_account=PAS0471                                     # If the scripts is submitted with another project, this will be updated (line below)
[[ -n $SLURM_JOB_ACCOUNT ]] && osc_account=$(echo "$SLURM_JOB_ACCOUNT" | tr "[:lower:]" "[:upper:]")
conda_path=/fs/project/PAS0471/jelmer/conda/nextflow
workflow_dir_base=workflows/nfcore-ampliseq             # Dir to download the workflow files to
container_dir=/fs/scratch/"$osc_account"/containers     # The workflow will download containers to this dir
work_dir=/fs/scratch/"$osc_account"/$USER/nfc-ampliseq  # 'work dir' for initial outputs (selected, final outputs go to the outdir)
profile="singularity"
resume=true && resume_opt="-resume"
version_only=false                                      # When true, just print tool & script version info and exit

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo "ABOUT:"
    echo "  - The workflow definition files will be downloaded by the script if not already present for the correct version"
    echo "  - All workflow parameter defaults in this script are the same as in the https://nf-co.re/ampliseq pipeline,"
    echo "    except that '--ignore_failed_trimming' is always used so the pipeline will keep running when some samples have too small filesizes after trimming."
    echo "  - Different from the Nextflow default, this script will try to 'resume' (rather than restart) a previous incomplete run by default"
    echo "  - This workflow can be used for both 16S and ITS data (use the '--its' flag for the latter)"
    echo "  - Not all nf-core/ampliseq parameters are present as options to this script, use '--more_opts' for parameters that aren't listed"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT \\"
    echo "    --more-args \"--extension '/*_R{1,2}.fastq.gz' --illumina_novaseq\""
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --input             <dir>   Input dir with FASTQ files"
    echo "  --outdir            <dir>   Output directory (will be created)"
    echo "  --FW_primer         <str>   Forward primer sequence"
    echo "  --RV_primer         <str>   Reverse primer sequence"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --its                       Use this flag if the data is ITS-derived (instead of 16S), in which case:"
    echo "                                - The default taxonomy will be changed to $ITS_taxonomy"
    echo "                                - The '--illumina_pe_its' and '--addsh' flags to the workflow will be used"
    echo "  --workflow_version  <str>   Nf-core ampliseq workflow version to use                    [default: $workflow_version]"
    echo "  --extension         <dir>   FASTQ file extension glob pattern                           [default: '/*_R{1,2}_001.fastq.gz']"
    echo "  --trunclenf         <int>   Truncate forward reads at this length                       [default: none => uses quality instead]"
    echo "  --trunclenr         <int>   Truncate reverse reads at this length                       [default: none => uses quality instead]"
    echo "  --min_len_asv       <int>   Minimum ASV length in basepairs                             [default: none => no length filtering]"
    echo "  --max_len_asv       <int>   Maximum ASV length in basepairs                             [default: none => no length filtering]"
    echo "  --sample_inference  <str>   Dada sample inference method                                [default: 'independent']"
    echo "                                - 'independent', 'pooled', or 'pseudo'"
    echo "  --dada_ref_taxonomy <str>   Reference taxonomy for dada                                 [16S default: Ampliseq's default ('silva=138' as of writing)]"
    echo "                                                                                          [ITS default: $ITS_taxonomy]"
    echo "  --filter_ssu        <str>   Quoted, comma-separated list of kingdoms to keep after Barrnap classification [default: no SSU filtering]"
    echo "                                - Recommended for 16S: 'bac', see https://nf-co.re/ampliseq/parameters#filter_ssu"
    echo "                                - Don't use this option for an ITS dataset!"
    echo "  --exclude_taxa      <str>   Quoted, comma-separated list of Kingdoms/taxa to remove     [default: 'mitochondria,chloroplast']"
    echo "                                after ASV tax. classification with DADA2 or QIIME2"
    echo "  --min_frequency     <int>   ASV must be present at least x times                        [default: 1]"
    echo "  --metadata          <file>  Metadata TSV file                                           [default: no metadata => no by-group stats]"
    echo "                                - At a minimum should have a column named 'id'"
    echo "                                - See https://docs.qiime2.org/2022.11/tutorials/metadata/"
    echo "  --metadata_category <str>   Metadata category/ies for downstream analysis               [default: none]"
    echo "  --metadata_category_barplot <str> Metadata category/ies for barplots                    [default: none / 'metadata_category' if that is specified]"
    echo "  --more_opts         <str>   Additional arguments, all in a single quoted string, to pass to the workflow"
    echo
    echo "NEXTFLOW AND UTILITY OPTIONS:"
    echo "  --restart                   Don't attempt to resume workflow run, but start over        [default: resume workflow]"
    echo "  --workflow_dir      <dir>   Top-level dir (no version number) with/for workflow files   [default: $workflow_dir_base]"
    echo "                                - If the correct version of the workflow is already present in this dir, it won't be downloaded again"
    echo "  --container_dir     <dir>   Directory with container images                             [default: $container_dir]"
    echo "                                - Required images will be downloaded here when not already present here" 
    echo "  --config            <file>  Additional config file                                      [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --profile            <str>  'Profile' to use from one of the config files               [default: $profile]"
    echo "  --work_dir           <dir>  Scratch (work) dir for the workflow                         [default: $work_dir]"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
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
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
fastq_dir=
outdir=
FW_primer=
RV_primer=
exclude_taxa= && exclude_taxa_opt=
extension= && extension_opt=
trunclenf= && trunclenf_opt=
trunclenr= && trunclenr_opt=
sample_inference= && sample_inference_opt=
dada_ref_taxonomy= && dada_ref_taxonomy_opt=
min_len_asv= && min_len_asv_opt=
max_len_asv= && max_len_asv_opt=
filter_ssu= && filter_ssu_opt=
config_file=
metadata= && metadata_opt=
metadata_category=
metadata_category_barplot=
min_frequency=1
threads=
more_opts=

# Parse command-line args
all_opts="$*"
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
        --workflow_version )            shift && workflow_version=$1 ;;
        --workflow_dir )                shift && workflow_dir_base=$1 ;;
        --container_dir )               shift && container_dir=$1 ;;
        --more_opts )                   shift && more_opts=$1 ;;
        --config | -config )            shift && config_file=$1 ;;
        --profile | -profile )          shift && profile=$1 ;;
        --work_dir | -work-dir )        shift && work_dir=$1 ;;
        --restart | -restart )          resume=false && resume_opt= ;;
        -h | --help )                   script_help; exit 0 ;;
        -v )                            script_version; exit 0 ;;
        --version )                     version_only=true ;;
        * )                             die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# Check arguments
[[ -z "$fastq_dir" ]] && die "No FASTQ dir specified, do so with -i/--input" "$all_opts"
[[ -z "$FW_primer" ]] && die "No forward primer specified, do so with --FW_primer" "$all_opts"
[[ -z "$RV_primer" ]] && die "No reverse primer specified, do so with --RV_primer" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -d "$fastq_dir" ]] && die "Input dir $fastq_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load software and set nr of threads
load_env "$conda_path"
nextflow_setup
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# ITS options
if [[ "$is_ITS" = true ]]; then
    ITS_opt="$ITS_opt_default"
    [[ -z "$dada_ref_taxonomy" ]] && dada_ref_taxonomy="$ITS_taxonomy"
fi

# Metadata options
if [[ -n "$metadata" ]]; then
    if [[ -n "$metadata_category" ]]; then
        [[ -z "$metadata_category_barplot" ]] && metadata_category_barplot="$metadata_category"
        metadata_opt="--metadata $metadata --metadata_category $metadata_category --metadata_category_barplot $metadata_category_barplot"
    else
        metadata_opt="--metadata $metadata"
    fi
fi

# Other options - only pass these args if non-default parameters were passed
[[ -n "$sample_inference" ]] && sample_inference_opt="--sample_inference $sample_inference"
[[ -n "$dada_ref_taxonomy" ]] && dada_ref_taxonomy_opt="--dada_ref_taxonomy $dada_ref_taxonomy"
[[ -n "$min_len_asv" ]] && min_len_asv_opt="--min_len_asv $min_len_asv"
[[ -n "$max_len_asv" ]] && max_len_asv_opt="--max_len_asv $max_len_asv"
[[ -n "$trunclenf" ]] && trunclenf_opt="--trunclenf $trunclenf"
[[ -n "$trunclenr" ]] && trunclenr_opt="--trunclenr $trunclenr"
[[ -n "$min_frequency" ]] && min_frequency_opt="--min_frequency $min_frequency"
[[ -n "$filter_ssu" ]] && filter_ssu_opt="--filter_ssu $filter_ssu"
[[ -n "$extension" ]] && extension_opt="--extension $extension"
[[ -n "$exclude_taxa" ]] && exclude_taxa_opt="--exclude_taxa $exclude_taxa"

# Build the config argument
OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_opt="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_opt="$config_opt -c ${config_file/,/ -c }"

# Other dirs
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
trace_dir="$outdir"/pipeline_info
workflow_dir_full="$workflow_dir_base"/${workflow_version//./_}

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:                 $all_opts"
echo
echo "INPUT AND OUTPUT:"
echo "Input FASTQ dir:                              $fastq_dir"
echo "Output dir:                                   $outdir"
echo "File extension:                               $extension"
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
[[ -n "$metadata_opt" ]] && echo "Metadata argument:                            $metadata_opt"
[[ -n "$more_opts" ]] && echo "Additional arguments:                         $more_opts"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run (if any):                 $resume"
echo "Container dir:                                $container_dir"
echo "Scratch (work) dir:                           $work_dir"
echo "Nextflow workflow dir:                        $workflow_dir_full"
echo "Config 'profile':                             $profile"
echo "Config file argument:                         $config_opt"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
set_threads "$IS_SLURM"
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

# Modify the config file so it has the correct OSC project/account
if [[ "$osc_account" != "PAS0471" ]]; then
    sed -i "s/--account=PAS0471/--account=$osc_account/" "$OSC_CONFIG"
fi

# Download workflow, if needed
if [[ ! -d "$workflow_dir_full" ]]; then
    mkdir -p "$(dirname "$workflow_dir_base")"
    echo -e "\n# Downloading workflow to $workflow_dir_base"
    nf-core download "$WORKFLOW_NAME" \
        --revision "$workflow_version" \
        --compress none \
        --container-system singularity \
        --container-cache-utilisation amend \
        --parallel-downloads "$threads" \
        --outdir "$workflow_dir_base" \
        --force
    echo
fi

# Run the tool
log_time "Starting the workflow.."
runstats $TOOL_BINARY \
    "$workflow_dir_full" \
    --input_folder "$fastq_dir" \
    $extension_opt \
    --outdir "$outdir" \
    --FW_primer "$FW_primer" \
    --RV_primer "$RV_primer" \
    --ignore_failed_trimming \
    $ITS_opt \
    $dada_ref_taxonomy_opt \
    $trunclenf_opt \
    $trunclenr_opt \
    $sample_inference_opt \
    $min_len_asv_opt \
    $max_len_asv_opt \
    $min_frequency_opt \
    $filter_ssu_opt \
    $exclude_taxa_opt \
    $metadata_opt \
    -work-dir "$work_dir" \
    -with-report "$trace_dir"/report.html \
    -with-trace "$trace_dir"/trace.txt \
    -with-timeline "$trace_dir"/timeline.html \
    -with-dag "$trace_dir"/dag.png \
    -ansi-log false \
    -profile "$profile" \
    $config_opt \
    $resume_opt \
    $more_opts

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
