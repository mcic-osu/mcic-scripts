#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=nfcore_rnaseq
#SBATCH --output=slurm-nfcore_rnaseq-%j.out


# FUNCTIONS --------------------------------------------------------------------
## Help function
Print_help() {
    echo
    echo "======================================================================================="
    echo "                             $0"
    echo "        Run the Nextflow-core RNAseq pipeline from https://nf-co.re/rnaseq"
    echo "======================================================================================="
    echo
    echo "USAGE:"
    echo "    sbatch $0 -i <samplesheet> -f <reference-fasta> -g <reference-annot> -o <outdir> [...]"
    echo "        (Always submit the script to the queue with 'sbatch' when wanting to start a workflow run)"
    echo "    bash $0 -h"
    echo "        (Run the script using 'bash' and the -h option to get help)"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i DIR      Sample sheet containing paths to FASTQ files and sample info"
    echo "                (See https://nf-co.re/rnaseq/3.9/usage#samplesheet-input)"
    echo "    -f FILE     Reference genome FASTA file"
    echo "    -g FILE     Reference genome annotation file (GTF/GFF - GTF preferred)"
    echo "    -o DIR      Output directory (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -r          Don't attempt to resume workflow run, but start over      [default: resume workflow]"
    echo "    -a STRING   Additional arguments to pass to 'nextflow run'"
    echo "    -c FILE     Additional config file                                    [default: none]"
    echo "                    - Settings in this file will override default settings"
    echo "                    - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                      (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "    -d DIR      Dir with the assembly workflow repository                 [default: 'workflows/nf-core-rnaseq']"
    echo "                    - If this dir isn't present, the workflow will be automatically downloaded."
    echo
    echo "OPTIONS YOU PROBABLY DON'T NEED TO USE:"
    echo "    -p STRING   Profile to use from one of the config files               [default: 'singularity']"
    echo "    -s DIR      Singularity container dir                                 [default: '/fs/project/PAS0471/containers']"
    echo "    -w DIR      Scratch (work) dir for the workflow                       [default: '/fs/scratch/PAS0471/\$USER/nfc_rnaseq']"
    echo "                    - This is where the workflow results will be stored"
    echo "                      before final results are copied to the output dir."
    echo
    echo "UTLITY OPTIONS:"
    echo "    -x          Turn on debugging/dry-run mode: print run information, but don't run"
    echo "    -h          Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "    sbatch $0 -i data/meta/samplesheet.csv -f data/ref/my.fa -a data/ref/my.gtf -o results/nfc_rnaseq"
    echo
    echo "OUTPUT:"
    echo "    - HTML file with summary of results: <outdir>/multiqc/star_salmon/multiqc_report.html"
    echo "    - Gene counts for use with DESeq2 in <outdir>/star_salmon/salmon.merged.gene_counts_length_scaled.rds"
    echo "      Use the script 'mcic-scripts/R_templates/nfcore-rnaseq_load-counts.R' to get started with that."
    echo
    echo "NFCORE RNASEQ WORKFLOW DOCUMENTATION:"
    echo " - https://nf-co.re/rnaseq "
    echo
}

Load_software() {
    module load python/3.6-conda5.2
    source activate /fs/project/PAS0471/jelmer/conda/nextflow

    ## Singularity container dir - any downloaded containers will be stored here;
    ## if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"

    ## Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting...\n" >&2
    exit 1
}


# CONTSTANTS AND DEFAULTS ------------------------------------------------------
## URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

## Option defaults
workflow_dir=workflows/nf-core-rnaseq
container_dir=/fs/project/PAS0471/containers
scratch_dir=/fs/scratch/PAS0471/$USER/nfc_rnaseq
profile=singularity
resume=true
debug=false


# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
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
        h) Print_help && exit 0 ;;
        \?) Print_help; Die "Invalid option $OPTARG" ;;
        :) Print_help; Die "Option -$OPTARG requires an argument" ;;
    esac
done


# SETUP ---------------------------------------------------------------
## Load software
Load_software

## Bash strict settings
set -ueo pipefail

## Check input
[[ "$samplesheet" = "" ]] && Die "Please specify a samplesheet with -i"
[[ "$ref_fasta" = "" ]] && Die "Please specify a genome FASTA file with -f"
[[ "$ref_annot" = "" ]] && Die "Please specify an annotation (GFF/GTF) file with -g"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f "$samplesheet" ]] && Die "Samplesheet $samplesheet does not exist"
[[ ! -f "$ref_fasta" ]] && Die "Reference FASTA file $ref_fasta does not exist"
[[ ! -f "$ref_annot" ]] && Die "Reference annotation file $ref_annot does not exist"

## Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    wget -q "$OSC_CONFIG_URL"
    osc_config=$(basename "$OSC_CONFIG_URL")
fi

## Build the config argument
config_arg="-c $osc_config"
if [[ "$config_file" != "" ]]; then
    config_arg="$config_arg -c ${config_file/,/ -c }"
fi

## Setup Nextflow arguments: resume option
if [[ "$resume" = true ]]; then
    resume_arg="-resume"
else
    resume_arg=""
fi

## Setup Nextflow arguments: annottation filetype
if [[ "$ref_annot" =~ .*\.gff3? ]]; then
    annot_arg="--gff $ref_annot"
elif [[ "$ref_annot" =~ .*\.gtf ]]; then
    annot_arg="--gtf $ref_annot"
else
    Die "Unknown annotation file format" >&2 && exit 1
fi

## Report
echo
echo "=========================================================================="
echo "              STARTING SCRIPT NCFORE_RNASEQ.SH"
date
echo "=========================================================================="
echo "Sample sheet:                      $samplesheet"
echo "Reference genome FASTA file:       $ref_fasta"
echo "Reference genome annotation file:  $ref_annot"
echo "Output dir:                        $outdir"
echo "Resume previous run:               $resume"
echo
echo "Container dir:                     $container_dir"
echo "Scratch (work) dir:                $scratch_dir"
echo "Dir with workflow files:           $workflow_dir"
echo "Config 'profile':                  $profile"
echo "Config file argument:              $config_arg"
[[ "$config_file" != "" ]] && echo "Additional config file:            $config_file"
[[ "$more_args" != "" ]] && echo "Additional arguments:              $more_args"
echo "=========================================================================="


# MAIN -------------------------------------------------------------------------
if [[ "$debug" = false ]]; then
    ## Make necessary dirs
    trace_dir="$outdir"/pipeline_info
    mkdir -p "$scratch_dir" "$container_dir" "$outdir" "$trace_dir"

    ## Remove old trace files
    [[ -f "$trace_dir"/report.html ]] && rm "$trace_dir"/report.html
    [[ -f "$trace_dir"/trace.txt ]] && rm "$trace_dir"/trace.txt
    [[ -f "$trace_dir"/timeline.html ]] && rm "$trace_dir"/timeline.html
    [[ -f "$trace_dir"/dag.png ]] && rm "$trace_dir"/dag.png

    ## Download workflow, if needed
    if [[ ! -d "$workflow_dir" ]]; then
        mkdir -p "$(dirname "$workflow_dir")"
        nf-core download rnaseq \
            --revision 3.9 \
            --compress none \
            --container singularity \
            --outdir "$workflow_dir"
    fi
fi

## Define the workflow command
run_workflow() {
    nextflow run "$workflow_dir"/workflow \
        --input "$samplesheet" \
        --outdir "$outdir" \
        --fasta "$ref_fasta" \
        $annot_arg \
        --aligner star_salmon \
        --remove_ribo_rna \
        --save_reference \
        --save_non_ribo_reads \
        --save_merged_fastq \
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
    echo -e "\n## Final command:"
    type run_workflow | sed '1,3d;$d'
    nextflow config "$workflow_dir" -profile "$profile"
else
    echo -e "\n## Starting the nextflow run..."
    eval run_workflow
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
echo
