#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=nf_dual
#SBATCH --output=slurm-nf_dual-%j.out


# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run the nf-core dual RNAseq pipeline."
  echo
  echo "Required options:"
  echo "    -i DIR      Input directory with FASTQ files"
  echo "    -f FILE     Input host FASTA file (NOTE: extension should be .fa or .fasta)"
  echo "    -F FILE     Input pathogen FASTA file (NOTE: extension should be .fa or .fasta)"
  echo "    -g FILE     Input host GFF file"
  echo "    -G FILE     Input pathogen GFF file"
  echo "    -o DIR      Output directory for workflow results"
  echo "    -c FILE     Config file with settings and SLURM execution details"
  echo "                NOTE: Reference FASTA and GFF files should be specified in this config file"
  echo
  echo "Other options:"
  echo "    -a STRING   Host GFF gene attribute/identifier     [default: 'gene_id]'"
  echo "    -A STRING   Pathogen GFF gene attribute/identifier [default: 'gene_id']"
  echo "    -p STRING   Profile                                [default: 'singularity']"
  echo "    -s DIR      Singularity container dir              [default: '/fs/project/PAS0471/containers']"
  echo "    -w DIR      Scratch (work) dir for the workflow    [default: '/fs/scratch/PAS0471/$USER/nf_dualrnaseq']"
  echo "    -r          Resume incomplete run                  [default: don't resume, i.e. start anew]"
  echo "    -x STRING   Sample ID glob to select only some samples"
  echo "    -h          Print this help message and exit"
  echo
  echo "Example command:"
  echo "    $0 -o results/nf_dual -p profile"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
container_dir=/fs/project/PAS0471/containers
scratch_dir=/fs/scratch/PAS0471/$USER/nf_dualrnaseq
profile=singularity

resume=false
resume_arg=""
fq_dir=""
outdir=""
config_org=""
sample_glob=""

fasta_host=""
fasta_pathogen=""
gff_host=""
gff_pathogen=""

attribute_host="gene_id"
attribute_pathogen="gene_id"

## Parse command-line options
while getopts 'i:o:f:F:g:G:c:x:p:w:s:a:A:rh' flag; do
  case "${flag}" in
  i) fq_dir="$OPTARG" ;;
  f) fasta_host="$OPTARG" ;;
  F) fasta_pathogen="$OPTARG" ;;
  g) gff_host="$OPTARG" ;;
  G) gff_pathogen="$OPTARG" ;;
  a) attribute_host="$OPTARG" ;;
  A) attribute_pathogen="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  c) config_org="$OPTARG" ;;
  w) scratch_dir="$OPTARG" ;;
  s) container_dir="$OPTARG" ;;
  p) profile="$OPTARG" ;;
  x) sample_glob="$OPTARG" ;;
  r) resume=true ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SOFTWARE SETUP ---------------------------------------------------------------
## Load Conda environment
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/nextflow

## Singularity container dir
export NXF_SINGULARITY_CACHEDIR="$container_dir"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

## Limit memory for Nextflow main process -https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
export NXF_OPTS='-Xms1g -Xmx7g'


# OTHER SETUP ------------------------------------------------------------------
## Bash strict settings
set -ueo pipefail

## Check input
[[ "$fq_dir" = "" ]] && echo "## ERROR: Please specify an input FASTQ dir with -i" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" && exit 1
[[ "$scratch_dir" = "" ]] && echo "## ERROR: Please specify a scratch (work) dir with -w" && exit 1
[[ "$config_org" = "" ]] && echo "## ERROR: Please specify a config file with -c" && exit 1
[[ "$fasta_host" = "" ]] && echo "## ERROR: Please specify a host FASTA file with -f" && exit 1
[[ "$fasta_pathogen" = "" ]] && echo "## ERROR: Please specify a pathogen FASTA file with -F" && exit 1
[[ "$gff_host" = "" ]] && echo "## ERROR: Please specify a host GFF file with -g" && exit 1
[[ "$gff_pathogen" = "" ]] && echo "## ERROR: Please specify a pathogen GFF file with -G" && exit 1

[[ ! -f "$fasta_host" ]] && echo "## ERROR: Input host FASTA file $fasta_host does not exist" && exit 1
[[ ! -f "$fasta_pathogen" ]] && echo "## ERROR: Input pathogen FASTA file $fasta_pathogen does not exist" && exit 1
[[ ! -f "$gff_host" ]] && echo "## ERROR: Input host GFF file $gff_host does not exist" && exit 1
[[ ! -f "$gff_pathogen" ]] && echo "## ERROR: Input pathogen GFF file $gff_pathogen does not exist" && exit 1
[[ ! -d "$fq_dir" ]] && echo "## ERROR: Input FASTQ dir $fq_dir does not exist" && exit 1
[[ ! -f "$config_org" ]] && echo "## ERROR: Input config file $config_org does not exist" && exit 1

## Process input options
[[ $resume = true ]] && resume_arg="-resume"

## Define files created by this script prior to running the workflow
config_final="$outdir"/config/$(basename "$config_org")

## Make output dirs
mkdir -p "$outdir"/config "$scratch_dir"

## Report
echo
echo "## Starting script nf_dual.sh"
date
echo
echo "## Input files:"
echo "## FASTQ dir:                       $fq_dir"
echo "## Host FASTA file:                 $fasta_host"
echo "## Pathogen FASTA file:             $fasta_pathogen"
echo "## Host GFF file:                   $gff_host"
echo "## Pathogen GFF file:               $gff_pathogen"
echo "## Input config file:               $config_org"
echo
echo "## Settings and output files:"
echo "## Host GFF attribute:              $attribute_host"
echo "## Pathogen GFF attribute:          $attribute_pathogen"
echo "## Output dir:                      $outdir"
echo "## Scratch dir:                     $scratch_dir"
echo "## Final config file:               $config_final"
echo "## Profile:                         $profile"
[[ $sample_glob != "" ]] && echo "## Sample glob:                     $sample_glob"
echo "## Resume run:                      $resume"
echo -e "-------------------------\n"


# MODIFY CONFIG ---------------------------------------------------------------
## Modify the config file: add paths to GFF files
#? Somehow GFF files need to be specified in the config file and not on the command-line -- the latter returns errors
sed -E "s@gff_host = .*@gff_host = \"$gff_host\"@" "$config_org" |
    sed -E "s@gff_pathogen = .*@gff_pathogen = \"$gff_pathogen\"@" \
    > "$config_final"

echo "## Showing the contents of the final config file:"
cat "$config_final"


# RUN THE NEXTFLOW WORKFLOW ----------------------------------------------------
echo -e "\n-------------------------------"
echo "## Starting the nextflow run..."
nextflow run \
    nf-core/dualrnaseq $resume_arg \
    -c "$config_final" \
    -w "$scratch_dir" \
    -profile "$profile" \
    --input "$fq_dir/$sample_glob*_R{1,2}*fastq.gz" \
    --fasta_host "$fasta_host" \
    --fasta_pathogen "$fasta_pathogen" \
    --gene_feature_gff_to_quantify_host "exon" \
    --host_gff_attribute "$attribute_host" \
    --gene_feature_gff_to_quantify_pathogen "exon" \
    --pathogen_gff_attribute "$attribute_pathogen" \
    --outdir "$outdir" \
    --genomes_ignore

#? Adding '--genomes_ignore' seems necessary to run the workflow but doesn't seem to do anything
#? For the workflow's own config file, see see ~/.nextflow/assets/nf-core/dualrnaseq/nextflow.config


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script nf_dual.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
