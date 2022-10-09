#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=nfc_rnaseq
#SBATCH --output=slurm-nfc_rnaseq-%j.out


# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run the Nextflow-core RNAseq pipeline from https://nf-co.re/rnaseq"
    echo
    echo "Required options:"
    echo "    -i DIR      Sample sheet containing paths to FASTQ files and sample info"
    echo "    -o DIR      Output directory for workflow results (will be created if needed)"
    echo
    echo "Other options:"
    echo "    -p STRING   Profile to use from one of the config files             [default: 'singularity']"
    echo "    -c FILE     (Additional) config file                                [default: none]"
    echo "                    - Any settings in this file will override the settings in the workflow's 'nextflow.config' file"
    echo "    -r          Don't attempt to resume workflow run, but start over    [default: resume]"
    echo "    -d DIR      Dir with the assembly workflow repository               [default: 'workflows/nf-core-rnaseq']"
    echo "                    - The repository can be downloaded using 'nf-core download rnaseq -r 3.9 -x singularity -c none --outdir nf-core-rnaseq'"
    echo "    -s DIR      Singularity container dir                               [default: '/fs/project/PAS0471/containers']"
    echo "                    - This is where any containers used in the workflow will be downloaded to"
    echo "    -w DIR      Scratch (work) dir for the workflow                     [default: '/fs/scratch/PAS0471/$USER/nfc_rnaseq']"
    echo "                    - This is where the workflow results will be stored before final results are copied to the specified output dir"
    echo "                    - This should preferable be a dir in OSC's scratch dir rather than in the main project dir"
    echo "    -a STRING   Additional arguments to pass to 'nextflow run'"
    echo "    -h          Print this help message and exit"
    echo
    echo "Example command:"
    echo "    sbatch $0 -i data/meta/samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/assembly"
}

## Option defaults
workflow_dir=workflows/nf-core-rnaseq/workflow
container_dir=/fs/project/PAS0471/containers
scratch_dir=/fs/scratch/PAS0471/$USER/nfc_rnaseq
profile=singularity
resume=true

debug=false #TODO add option
samplesheet=""
ref_annot=""
ref_fasta=""
outdir=""
config_file=""
more_args=""

## Parse command-line options
while getopts 'i:g:f:o:c:d:w:s:p:a:hrx' flag; do
    case "${flag}" in
        i) samplesheet="$OPTARG" ;;
        g) ref_annot="$OPTARG" ;;
        f) ref_fasta="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) workflow_dir="$OPTARG" ;;
        w) scratch_dir="$OPTARG" ;;
        s) container_dir="$OPTARG" ;;
        p) profile="$OPTARG" ;;
        c) config_file="$OPTARG" ;;
        r) resume=false ;;
        a) more_args="$OPTARG" ;;
        x) debug=true ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SOFTWARE SETUP ---------------------------------------------------------------
## Load Conda environment
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/nextflow

## Get the OSC config file
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
OSC_CONFIG_FILE="$SCRIPTPATH"/../nextflow/osc.config
if [[ ! -f "$OSC_CONFIG_FILE" ]]; then
    echo "## Cloning mcic-scripts repo..."
    [[ ! -d "mcic-scripts" ]] && git clone https://github.com/mcic-osu/mcic-scripts.git
    OSC_CONFIG_FILE=mcic-scripts/nextflow/osc.config
fi
[[ ! -f "$OSC_CONFIG_FILE" ]] && echo "## No OSC config file" >&2 && exit 1
config_arg="-c $OSC_CONFIG_FILE"

## Singularity container dir - any downloaded containers will be stored here;
## if the required container is already there, it won't be re-downloaded
export NXF_SINGULARITY_CACHEDIR="$container_dir"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

## Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
export NXF_OPTS='-Xms1g -Xmx7g'


# OTHER SETUP ------------------------------------------------------------------
## Bash strict settings
set -ueo pipefail

## Check input
[[ "$samplesheet" = "" ]] && echo "## ERROR: Please specify a samplesheet with -i" && exit 1
[[ "$ref_fasta" = "" ]] && echo "## ERROR: Please specify a genome FASTA file with -f" && exit 1
[[ "$ref_annot" = "" ]] && echo "## ERROR: Please specify an annotation (GFF/GTF) file with -g" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" && exit 1
[[ ! -f "$samplesheet" ]] && echo "## ERROR: Samplesheet $samplesheet does not exist" && exit 1
[[ ! -f "$ref_fasta" ]] && echo "## ERROR: Reference FASTA file $ref_fasta does not exist" && exit 1
[[ ! -f "$ref_annot" ]] && echo "## ERROR: Reference annotation file $ref_annot does not exist" && exit 1

## Setup Nextflow arguments
#if [[ "$config_file" != "" ]]; then
#    config_arg="-c $config_file"
#else
#    config_arg=""
#fi

if [[ "$resume" = true ]]; then
    resume_arg="-resume"
else
    resume_arg=""
fi

if [[ "$ref_annot" =~ .*\.gff3? ]]; then
    annot_arg="--gff $ref_annot"
elif [[ "$ref_annot" =~ .*\.gtf ]]; then
    annot_arg="--gtf $ref_annot"
else
    echo "## ERROR: Unknown annotation file format" >&2 && exit 1
fi

## Report
echo
echo "## Starting script nfc_rnaseq.sh"
date
echo
echo "## Sample sheet:                      $samplesheet"
echo "## Reference genome FASTA file:       $ref_fasta"
echo "## Reference genome annotation file:  $ref_annot"
echo "## Output dir:                        $outdir"
echo
echo "## Container dir:                     $container_dir"
echo "## Scratch (work) dir:                $scratch_dir"
echo "## Dir with workflow files:           $workflow_dir"
echo "## Config 'profile':                  $profile"
echo "## Config file arg:                   $config_arg"
echo "## Resume previous run:               $resume"
[[ "$more_args" != "" ]] && echo "## Additional arguments:              $more_args"
[[ "$config_file" != "" ]] && echo "## Config file:                       $config_file"
echo -e "-------------------------\n"


# MAIN -------------------------------------------------------------------------
## Make necessary dirs
trace_dir="$outdir"/pipeline_info
mkdir -p "$scratch_dir" "$container_dir" "$outdir" "$trace_dir"

## Download workflow, if needed
if [[ ! -d "$workflow_dir" ]]; then
    mkdir -p "$(dirname "$workflow_dir")"
    #git clone "$WORKFLOW_URL" "$workflow_dir"
    nf-core download rnaseq -r 3.9 -x singularity -c none --outdir "$workflow_dir"
fi

## Define the workflow command
run_workflow() {
    nextflow run "$workflow_dir" \
        --input "$samplesheet" \
        --outdir "$outdir" \
        --fasta "$ref_fasta" \
        $annot_arg \
        --save_reference \
        --aligner star_salmon \
        --remove_ribo_rna \
        --save_non_ribo_reads \
        -work-dir "$scratch_dir" \
        -with-report "$trace_dir"/report.html \
        -with-trace "$trace_dir"/trace.txt \
        -with-timeline "$trace_dir"/timeline.html \
        -with-dag "$trace_dir"/dag.png \
        -ansi-log false \
        -profile "$profile" $config_arg $resume_arg $more_args
}

## Run/echo the workflow
if [[ "$debug" = true ]]; then
    echo "## Final command:"
    type run_workflow | sed '1,3d;$d'
    nextflow config "$workflow_dir" -profile "$profile"
fi

if [[ "$debug" != true ]]; then
    echo -e "\n## Starting the nextflow run..."
    eval run_workflow
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script nfc_rnaseq.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
