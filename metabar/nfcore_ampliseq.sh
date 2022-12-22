#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=nfcore_ampliseq
#SBATCH --output=slurm-nfcore_ampliseq-%j.out

# ==============================================================================
#                                FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================================="
    echo "                             $0"
    echo "        Run the Nextflow-core metabarcoding pipeline from https://nf-co.re/ampliseq"
    echo "======================================================================================="
    echo
    echo "ABOUT:"
    echo "  - All workflow parameter defaults in this script are the same as in the https://nf-co.re/ampliseq pipeline"
    echo "  - Different from the Nextflow default, this script will try to 'resume' (rather than restart) a previous incomplete run by default"
    echo "  - This workflow can be used for both 16S and ITS data (use the '--its' flag for the latter)"
    echo "  - Not all nf-core/ampliseq parameters are present as options to this script, use '--more_args' for parameters that aren't listed"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <indir> -o <outdir> --FW_primer <forward primer> --RV_primer <reverse primer> [...]"
    echo "    (Always submit the script to the queue with 'sbatch' when wanting to start a workflow run)"
    echo "  bash $0 -h"
    echo "    (Run the script using 'bash' and the -h option to get help)"
    echo
    echo "NF-CORE/AMPLISEQ WORKFLOW PARAMETERS:"
    echo "Note: These are all the same as the parameters specified at https://nf-co.re/ampliseq/parameters: see there for details"
    echo "Note: To pass nf-core/ampliseq parameters that are not listed below, use the 'more_args' option."
    echo "  --input             <dir>   Input FASTQ dir                         [REQUIRED]"
    echo "  --outdir            <dir>   Output directory (will be created)      [REQUIRED]"
    echo "  --FW_primer         <str>   Forward primer sequence                 [REQUIRED]"
    echo "  --RV_primer         <str>   Reverse primer sequence                 [REQUIRED]"
    echo "  --extension         <dir>   FASTQ file extension glob pattern       [default: '/*_R{1,2}_001.fastq.gz']"
    echo "  --trunclenf         <int>   Truncate forward reads at this length   [default: none => uses quality instead]"
    echo "  --trunclenr         <int>   Truncate reverse reads at this length   [default: none => uses quality instead]"
    echo "  --min_len_asv       <int>   Minimum ASV length in basepairs         [default: none => no length filtering]"
    echo "  --max_len_asv       <int>   Maximum ASV length in basepairs         [default: none => no length filtering]"
    echo "  --sample_inference  <str>   Dada sample inference method            [default: 'independent']"
    echo "                                - 'independent', 'pooled', or 'pseudo'"
    echo "  --dada_ref_taxonomy <str>   Reference taxonomy for dada             [default: 'silva=138' for 16S / 'unite-fungi=8.3' for ITS]"
    echo "  --filter_ssu        <str>   #TODO"
    echo "  --min_frequency     <int>   ASV must be present at least x times    [default: 1]"
    echo "  --metadata          <file>  Metadata TSV file                       [default: no metadata => no by-group stats]"
    echo "                                - At a minimum should have a column named 'id'"
    echo "                                - See https://docs.qiime2.org/2022.11/tutorials/metadata/"
    echo "  --metadata_category <str>   Metadata category/ies for downstream analysis [default: none]"
    echo "  --metadata_category_barplot <str> Metadata category/ies for barplots [default: none / 'metadata_category' if that is specified]"
    echo
    echo "OPTIONS OF THIS SHELL SCRIPT:"
    echo "  --its                       Use this flag if the data is ITS, in which case:"
    echo "                                - The default taxonomy will be changed to 'unite-fungi=8.3'"
    echo "                                - The '--illumina_pe_its' and '--addsh' flags to the workflow will be used"
    echo "  --more_args         <str>   Additional arguments to pass to 'nextflow run', pass as a single quoted string"
    echo "  --nf_file           <file>  Main .nf workflow definition file       [default: 'workflows/nfcore-ampliseq/workflow/main.nf']"
    echo "                                - If this file isn't present, the workflow will be automatically downloaded."
    echo "  --debug                     Turn on debugging mode: print all commands"
    echo "  --dryrun                    Dry run: print info and main commands but don't run"
    echo "  -h/--help                   Print this help message and exit"
    echo 
    echo "NEXTFLOW OPTIONS:"
    echo "  -r/-no-resume               Don't attempt to resume workflow run, but start over        [default: resume workflow]"
    echo "  -c/-config          <file>  Additional config file                                      [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  -profile            <str>   'Profile' to use from one of the config files               [default: 'singularity']"
    echo "  -work-dir           <dir>   Scratch (work) dir for the workflow                         [default: '/fs/scratch/PAS0471/\$USER/nfcore-ampliseq']"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT"
    echo "  sbatch $0 -i data/fastq -o results/ampliseq --FW_primer GTGTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT \\"
    echo "    --more-args \"--extension '/*_R{1,2}.fastq.gz' --illumina_novaseq\""
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
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
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
# Option defaults - workflow parameters
extension='/*_R{1,2}_001.fastq.gz'     # Same as nfcore/ampliseq default
is_ITS=false && ITS_arg=""

# Option defaults - Nextflow
nf_file=workflows/nfcore-ampliseq/workflow/main.nf
container_dir=/fs/project/PAS0471/containers
work_dir=/fs/scratch/PAS0471/$USER/nfcore-ampliseq
profile="singularity"
resume=true && resume_arg="-resume"

# Option defaults - utitilites
debug=false
dryrun=false && e=""
slurm=true

# URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here


# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
# Placeholder defaults
fastq_dir=""
outdir=""
FW_primer=""
RV_primer=""
trunclenf=""
trunclenr=""
sample_inference=""
dada_ref_taxonomy=""
min_len_asv="" && min_len_asv_arg=""
max_len_asv="" && max_len_asv_arg=""
min_frequency=1
filter_ssu=""
config_file=""
metadata="" && metadata_arg=""
metadata_category=""
metadata_category_barplot=""
more_args=""
sample_inference="pseudo"

# Parse command-line options
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
        --metadata )                    shift && metadata=$1 ;;
        --metadata_category )           shift && metadata_category=$1 ;;
        --metadata_category_barplot )   shift && metadata_category_barplot=$1 ;;
        --nf_file )                     shift && nf_file=$1 ;;
        --container_dir )               shift && container_dir=$1 ;;
        --more_args )                   shift && more_args=$1 ;;
        -c | -config )                  shift && config_file=$1 ;;
        -p | -profile )                 shift && profile=$1 ;;
        -w | -work-dir )                shift && work_dir=$1 ;;
        -r | -no-resume )               resume=false ;;
        --debug )                       debug=true ;;
        --dryrun )                      dryrun=true ;;
        -h | --help )                   Print_help; exit ;;
        * )                             Print_help; Die "Invalid option $1";;
    esac
    shift
done


# ==============================================================================
#                              OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

# Load Conda environment
[[ "$dryrun" = false ]] && Load_software

# Bash strict settings
set -ueo pipefail

# Check input
[[ "$fastq_dir" = "" ]] && Die "Please specify an input FASTQ dir with -i/--input" "$all_args"
[[ "$FW_primer" = "" ]] && Die "Please specify a forward primer sequence with --FW_primer" "$all_args"
[[ "$RV_primer" = "" ]] && Die "Please specify a reverse primer sequence with --RV_primer" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o" "$all_args"

# ITS options
if [[ "$is_ITS" = true ]]; then
    ITS_arg="--illumina_pe_its --addsh"
    [[ "$dada_ref_taxonomy" = "" ]] && dada_ref_taxonomy='unite-fungi=8.3'
fi

# Metadata options
if [[ "$metadata" != "" ]]; then
    if [[ "$metadata_category" != "" ]]; then
        [[ "$metadata_category_barplot" = "" ]] && metadata_category_barplot="$metadata_category"
        metadata_arg="--metadata $metadata --metadata_category $metadata_category --metadata_category_barplot $metadata_category_barplot"
    else
        metadata_arg="--metadata $metadata"
    fi
fi

# Other options - only pass these args if non-default parameters were passed
[[ "$sample_inference" != "" ]] && sample_inference_arg="--sample_inference $sample_inference"
[[ "$dada_ref_taxonomy" != "" ]] && dada_ref_taxonomy_arg="--dada_ref_taxonomy $dada_ref_taxonomy"
[[ "$min_len_asv" != "" ]] && min_len_asv_arg="--min_len_asv $min_len_asv"
[[ "$max_len_asv" != "" ]] && max_len_asv_arg="--max_len_asv $max_len_asv"
[[ "$trunclenf" != "" ]] && trunclenf_arg="--trunclenf $trunclenf"
[[ "$trunclenr" != "" ]] && trunclenr_arg="--trunclenr $trunclenr"
[[ "$min_frequency" != "" ]] && min_frequency_arg="--min_frequency $min_frequency"
[[ "$filter_ssu" != "" ]] && filter_ssu_arg="--filter_ssu $filter_ssu"

# Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    wget -q "$OSC_CONFIG_URL"
    osc_config=$(basename "$OSC_CONFIG_URL")
fi

# Build the config argument
config_arg="-c $osc_config"
if [[ "$config_file" != "" ]]; then
    config_arg="$config_arg -c ${config_file/,/ -c }"
fi

# Setup Nextflow arguments: resume option
[[ "$resume" = false ]] && resume_arg=""

# Other output dirs
trace_dir="$outdir"/pipeline_info

# Report
echo
echo "=========================================================================="
echo "              STARTING SCRIPT NCFORE_AMPLISEQ.SH"
echo "=========================================================================="
date
echo "All arguments to this script:     $all_args"
echo "Input FASTQ dir:                  $fastq_dir"
echo "Output dir:                       $outdir"
echo "Forward primer:                   $FW_primer"
echo "Reverse primer:                   $RV_primer"
echo "Dada reference taxonomy:          $dada_ref_taxonomy"
[[ "$extension" != "" ]] && echo "File extension:                   $extension"
[[ "$sample_inference" != "" ]] && echo "Sample inference method:          $sample_inference"
[[ "$min_len_asv" != "" ]] && echo "Min. ASV length:                  $min_len_asv"
[[ "$max_len_asv" != "" ]] && echo "Max. ASV length:                  $max_len_asv"
[[ "$trunclenf" != "" ]] && echo "Truncate forward reads at:        $trunclenf"
[[ "$trunclenr" != "" ]] && echo "Truncate reverse reads at:        $trunclenr"
[[ "$filter_ssu" != "" ]] && echo "Filter SSU:                       $filter_ssu"
[[ "$min_frequency" != "" ]] && echo "Min. ASV frequency:               $min_frequency"
[[ "$metadata" != "" ]] && echo "Metadata file:                    $metadata"
[[ "$metadata_category" != "" ]] && echo "Metadata categories:              $metadata_category"
[[ "$metadata_category_barplot" != "" ]] && echo "Metadata categories for barplots: $metadata_category_barplot"
[[ "$metadata_arg" != "" ]] && echo "Metadata argument:                $metadata_arg"
[[ "$more_args" != "" ]] && echo "Additional arguments:             $more_args"
echo
echo "Resume previous run:              $resume"
echo "Container dir:                    $container_dir"
echo "Scratch (work) dir:               $work_dir"
echo "Nextflow workflow file:           $nf_file"
echo "Config 'profile':                 $profile"
echo "Config file argument:             $config_arg"
[[ "$config_file" != "" ]] && echo "Additional config file:           $config_file"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="
echo

# ==============================================================================
#                              RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    # Make necessary dirs
    mkdir -p "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

    # Remove old trace files
    [[ -f "$trace_dir"/report.html ]] && rm "$trace_dir"/report.html
    [[ -f "$trace_dir"/trace.txt ]] && rm "$trace_dir"/trace.txt
    [[ -f "$trace_dir"/timeline.html ]] && rm "$trace_dir"/timeline.html
    [[ -f "$trace_dir"/dag.png ]] && rm "$trace_dir"/dag.png

    # Download workflow, if needed
    if [[ ! -f "$nf_file" ]]; then
        workflow_dir="$(dirname "$(dirname "$nf_file")")"
        mkdir -p "$(dirname "$workflow_dir")"
        echo "# Downloading workflow to $workflow_dir"
        nf-core download ampliseq \
            --revision 2.4.1 \
            --compress none \
            --container singularity \
            --outdir "$workflow_dir"
        echo
    fi
fi

# Define the workflow command
echo -e "# Starting the workflow...\n"
[[ "$dryrun" = false ]] && set -o xtrace

${e}Time nextflow run \
    "$nf_file" \
    --input "$fastq_dir" \
    --extension "$extension" \
    --outdir "$outdir" \
    --FW_primer "$FW_primer" \
    --RV_primer "$RV_primer" \
    $ITS_arg \
    "${dada_ref_taxonomy_arg[@]}" \
    $trunclenf_arg \
    $trunclenr_arg \
    $sample_inference_arg \
    $min_len_asv_arg \
    $max_len_asv_arg \
    $min_frequency_arg \
    $filter_ssu_arg \
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

[[ "$debug" = false ]] && set +o xtrace


# ==============================================================================
#                              WRAP UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
