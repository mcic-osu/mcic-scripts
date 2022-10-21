#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=nf_resfinder
#SBATCH --output=slurm-nf_resfinder-%j.out

# FUNCTIONS --------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "                            $0: Run a ResFinder pipeline"
    echo "                      This pipeline will run AMR prediction with ResFinder"
    echo "=================================================================================================="
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i DIR        Input directory with FASTQ files"
    echo "    -s STRING     Species name (put in quotes, e.g. 'Salmonella enterica')"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "------------------"
    echo "    -o DIR      Output directory for workflow results                   [default: 'results/nf_resfinder']"
    echo "    -u          Don't infer point mutations                             [default: do infer point mutations]"
    echo "    -n DIR      Workflow definition file (a '*.nf' file)                [default: '/fs/project/PAS0471/jelmer/assist/2022-09_alejandra/workflows/resfinder/main.nf']"
    echo "    -g STRING   Sample ID glob to select only some samples/files        [default: '*R{1,2}*.fastq.gz']"
    echo "                    - Modify if your sample names don't adhere to the default, e.g. if they are 'fq.gz'"
    echo "                    - Can also use this to subset samples, e.g. 'sampleA*R{1,2}*.fastq.gz'"
    echo "    -r          Don't attempt to resume the workflow, but start over    [default: resume where you left off]"
    echo
    echo "OPTIONS YOU PROBABLY DON'T NEED TO USE:"
    echo "------------------"
    echo "    -t DIR      Singularity container dir                               [default: '/fs/project/PAS0471/containers']"
    echo "                    - This is where any containers used in the workflow will be downloaded to"
    echo "    -w DIR      Scratch (work) dir for the workflow                     [default: '/fs/scratch/PAS0471/$USER/mlst']"
    echo "                    - This is where the workflow results will be stored before final results are copied to the specified output dir"
    echo "                    - This should preferable be a dir in OSC's scratch dir rather than in the main project dir"
    echo "    -p STRING   Profile from any of the config files to use             [default: 'conda']"
    echo "    -c FILE     Additional config file(s)"
    echo "                    - Any settings in this file will override settings in default config files"
    echo "                    - Use a comma-separated list when supplying multiple files" 
    echo "                    - The mcic-scripts OSC config will always be used, https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config"
    echo "    -a STRING   Additional arguments to pass to 'nextflow run'"
    echo
    echo "UTILITY OPTIONS:"
    echo "------------------"
    echo "    -x          Turn on debugging/dry-run mode: print run information, but don't run commands"
    echo "    -h          Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "    sbatch $0 -i data/fastq"
    echo
    echo "DOCUMENTATION:"
    echo "------------------"
    echo "  - https://bitbucket.org/genomicepidemiology/resfinder/src/master/"
    echo
}

## Load the software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/nextflow
    
    ## Singularity container dir - any downloaded containers will be stored here
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    ## Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}


# OPTION DEFAULTS AND HARDCODED VALUES -----------------------------------------
## URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

## Option defaults
fastq_pattern='*R{1,2}*.fastq.gz'
outdir="results/nf_resfinder"
nextflow_file="/fs/project/PAS0471/jelmer/assist/2022-09_alejandra/workflows/resfinder/main.nf"
profile="conda"
container_dir=/fs/project/PAS0471/containers
scratch_dir=/fs/scratch/PAS0471/$USER/nf_resfinder
resume=true
debug=false
point_mut=true

# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
indir=""
species=""
config_file=""
more_args=""

## Parse command-line options
while getopts 'i:o:s:n:w:t:p:c:g:a:uxhr' flag; do
    case "${flag}" in
        i) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        s) species="$OPTARG" ;;
        n) nextflow_file="$OPTARG" ;;
        w) scratch_dir="$OPTARG" ;;
        t) container_dir="$OPTARG" ;;
        p) profile="$OPTARG" ;;
        c) config_file="$OPTARG" ;;
        g) fastq_pattern="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        u) point_mut=false ;;
        r) resume=false ;;
        x) debug=true ;;
        h) Help && exit 0 ;;
        \?) echo "$0: ERROR: Invalid option" >&2 && exit 1 ;;
        :) echo "$0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Load Conda environment
Load_software

## Bash strict settings
set -ueo pipefail

## Check input
[[ "$indir" = "" ]] && echo "$0: ERROR: Please specify an input dir with -i" && exit 1
[[ "$species" = "" ]] && echo "$0: ERROR: Please specify a species with -s" && exit 1
[[ ! -d "$indir" ]] && echo "$0: ERROR: Input dir $indir does not exist" && exit 1
[[ ! -f "$nextflow_file" ]] && echo "$0: ERROR: Nextflow file $nextflow_file does not exist" && exit 1

## Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    wget "$OSC_CONFIG_URL"
    osc_config=$(basename "$OSC_CONFIG_URL")
fi

## Build the config argument
config_arg="-c $osc_config"
if [[ "$config_file" != "" ]]; then
    config_arg="$config_arg -c ${config_file/,/ -c }"
fi

if [[ "$resume" = true ]]; then resume_arg="-resume"; else resume_arg=""; fi

if [[ "$point_mut" = true ]]; then
    point_arg="--resfinder_point_mutation"
else
    point_arg=""
fi

## Define trace output dir
trace_dir="$outdir"/pipeline_info

## Report
echo
echo "=========================================================================="
echo "                     STARTING SCRIPT NF_RESFINDER.SH"
date
echo "=========================================================================="
echo "## Input dir:                       $indir"
echo "## FASTQ file pattern:              $fastq_pattern"
echo "## Species:                         $species"
echo "## Output dir:                      $outdir"
echo "## Nextflow workflow file:          $nextflow_file"
echo
echo "## Container dir:                   $container_dir"
echo "## Scratch (work) dir:              $scratch_dir"
echo "## Config 'profile':                $profile"
echo "## Resume previous run:             $resume"
[[ "$more_args" != "" ]] && echo "## Additional arguments:            $more_args"
[[ "$config_file" != "" ]] && echo "## Config file:                     $config_file"
echo "=========================================================================="
echo


# MAIN -------------------------------------------------------------------------
if [[ "$debug" = false ]]; then
    ## Make necessary dirs
    mkdir -p "$scratch_dir" "$container_dir" "$outdir" "$trace_dir"

    ## Remove old trace files
    [[ -f "$trace_dir"/report.html ]] && rm "$trace_dir"/report.html
    [[ -f "$trace_dir"/trace.txt ]] && rm "$trace_dir"/trace.txt
    [[ -f "$trace_dir"/timeline.html ]] && rm "$trace_dir"/timeline.html
    [[ -f "$trace_dir"/dag.png ]] && rm "$trace_dir"/dag.png
    echo
fi

## Define the workflow run command
command="nextflow run $nextflow_file \
    --input_dir $indir \
    --fastq_pattern '$fastq_pattern' \
    --output_dir $outdir \
    --species '$species' \
    $point_arg \
    -work-dir $scratch_dir \
    -ansi-log false \
    -with-report $trace_dir/report.html \
    -with-trace $trace_dir/trace.txt \
    -with-timeline $trace_dir/timeline.html \
    -with-dag $trace_dir/dag.png \
    -profile $profile \
    $config_arg \
    $resume_arg \
    $more_args"

#? - Note that the docs say the adapter argument is 'adapter_sequences', but it is actually 'adapter_file'

## Run the workflow
echo "## Starting the workflow using the following command:"
echo
echo "$command" | tr -s " "
echo
if [[ "$debug" = false ]]; then
    eval $command
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
if [[ "$debug" = false ]]; then
    echo "## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
    echo
fi
echo "## Done with script"
date
