#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=bactopia
#SBATCH --output=slurm-bactopia-%j.out

# FUNCTIONS --------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "$0: Run the Bactopia pipeline"
    echo "=================================================================================================="
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i DIR      Sample sheet (FOFN) containing paths to FASTQ files"
    echo "    -o DIR      Output directory (will be created if needed)"
    echo "    -s STRING   Focal species, e.g. 'Salmonella enterica' (make sure to quote!)"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -r          Don't attempt to resume workflow run, but start over      [default: resume workflow]"
    echo "    -d DIR      Bactopia datasets directory                               [default: none"
    echo "    -p STRING   Profile to use from one of the config files               [default: 'singularity']"
    echo "    -c FILE     Additional config file                                    [default: mcic-scripts OSC config file 'osc.config']"
    echo "                    - osc.config is at https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config"
    echo "    -t DIR      Container dir                                             [default: /fs/project/PAS0471/containers]"
    echo "    -a STRING   Additional arguments to pass to Bactopia"
    echo
    echo "UTILITY OPTIONS:"
    echo "------------------"
    echo "    -x          Turn on debugging/dry-run mode: print run information, but don't run commands"
    echo "    -v          Print the version and exit"
    echo "    -h          Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "    sbatch $0 -i data/meta/samplesheet.tsv -s 'Salmonella enterica' -o results/bactopia"
    echo
    echo "OUTPUT:"
    echo "------------------"
    echo "      - ...."
    echo
    echo "BACTOPIA INFORMATION:"
    echo "------------------"
    echo " - Documentation: https://bactopia.github.io/ "
    echo " - Paper:         https://journals.asm.org/doi/10.1128/mSystems.00190-20"
}

## Print the version
Print_version() {
    load_software
    bactopia --version
}

## Load software
Load_software() {
    ## Load Conda environment
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate; done
    source activate /fs/project/PAS0471/jelmer/conda/bactopia

    ## Singularity container dir - any downloaded containers will be stored here;
    ## if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"

    ## Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}


# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
## Hardcoded variables
QUEUE_SIZE=100
MAX_CPUS=48
MAX_TIME=1440    # In hours
MAX_MEMORY=128   # In GB

## Option defaults
container_dir=/fs/project/PAS0471/containers
profile=singularity
resume=true
debug=false

datasets_dir=""
species=""
samplesheet=""
outdir=""
config_file=""
more_args=""

## Parse command-line options
while getopts 'i:o:s:d:c:p:t:a:rvhx' flag; do
    case "${flag}" in
        i) samplesheet="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        s) species="$OPTARG" ;;
        d) datasets_dir="$OPTARG" ;;
        p) profile="$OPTARG" ;;
        c) config_file="$OPTARG" ;;
        t) container_dir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        r) resume=false ;;
        x) debug=true ;;
        v) Print_version && exit 0 ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------
## Load software
Load_software

## Bash strict settings
set -ueo pipefail

## Check input
[[ "$samplesheet" = "" ]] && echo "## ERROR: Please specify a samplesheet with -i" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" && exit 1
[[ "$species" = "" ]] && echo "## ERROR: Please specify a species with -s" && exit 1
[[ ! -f "$samplesheet" ]] && echo "## ERROR: Samplesheet $samplesheet does not exist" && exit 1

## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT BACTOPIA.SH"
date
echo -e "==========================================================================\n"

## Get the OSC config file
if [[ "$config_file" = "" ]]; then
    echo "## Using the mcisc-scripts OSC config file"

    SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
    OSC_CONFIG_FILE="$PWD"/"$SCRIPTPATH"/../nextflow/osc.config
    
    if [[ ! -f "$OSC_CONFIG_FILE" ]]; then
        [[ ! -d "mcic-scripts" ]] && git clone https://github.com/mcic-osu/mcic-scripts.git
        OSC_CONFIG_FILE="$PWD"/mcic-scripts/nextflow/osc.config
        [[ ! -f "$OSC_CONFIG_FILE" ]] && echo "## ERROR: No OSC config file" >&2 && exit 1
    fi
    
    config_arg="--nfconfig $OSC_CONFIG_FILE"
fi

## Make paths absolute
[[ ! "$samplesheet" =~ ^/ ]] && samplesheet="$PWD"/"$samplesheet"

## Setup arguments: resume option
if [[ "$resume" = true ]]; then
    resume_arg="-resume"
else
    resume_arg=""
fi

## Setup arguments: datasets
if [[ "$datasets_dir" != "" ]]; then
    [[ ! "$datasets_dir" =~ ^/ ]] && datasets_dir="$PWD"/"$datasets_dir"
    datasets_arg="--datasets $datasets_dir"
else
    datasets_arg=""
fi

## Split outdir into 'outdir' and 'run_name'
outdir_base=$(dirname "$outdir")
run_name=$(basename "$outdir")

## Report
echo "## Sample sheet:                      $samplesheet"
echo "## Output dir:                        $outdir"
echo "## Species:                           $species"
[[ "$datasets_dir" != "" ]] && echo "## Bactopia datasets dir:             $datasets_dir"
echo
echo "## Config 'profile':                  $profile"
[[ "$config_file" != "" ]] && echo "## Additional config file:                       $config_file"
echo "## Resume previous run:               $resume"
[[ "$more_args" != "" ]] && echo "## Additional arguments:              $more_args"
echo
echo "Print the contents of the sample sheet:"
cat -n "$samplesheet"
echo -e "-------------------------\n"


# MAIN -------------------------------------------------------------------------
if [[ "$debug" = false ]]; then
    ## Make output dir
    mkdir -p "$outdir"/logs

    ## Move into output dir
    echo "## Changing into dir $outdir_base..."
    cd "$outdir_base" || exit 1
    echo

    ## Remove any previous trace files
    trace_file=$run_name/nf-reports/bactopia-trace.txt
    dag_file=$run_name/nf-reports/bactopia-dag.svg
    report_file=$run_name/nf-reports/bactopia-report.html
    timeline_file=$run_name/nf-reports/bactopia-timeline.html
    rm -vf "$timeline_file" "$report_file" "$dag_file" "$trace_file"
    echo
fi

## Define the workflow command
command="bactopia \
    --samples $samplesheet \
    --outdir ./ \
    --run_name $run_name \
    --force \
    -profile $profile \
    --species '$species' \
    --max_cpus $MAX_CPUS \
    --max_time $MAX_TIME \
    --max_memory $MAX_MEMORY \
    -qs $QUEUE_SIZE \
    $resume_arg \
    $datasets_arg \
    $config_arg \
    $more_args"

#? Apparently there is no setting to change the work dir
                    
## Report
echo "## Running Bactopia using the following command:"
echo "$command" | tr -s " "
echo

## Run Bactopia
[[ "$debug" = false ]] && eval $command


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
if [[ "$debug" = false ]]; then
    echo "## Ran Bactopia version:"
    Print_version | tee logs/version.txt
    echo "## Listing files in the output dir:"
    ls -lh
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
    echo
fi
echo "## Done with script bactopia.sh"
date
echo
