#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=trinity
#SBATCH --output=slurm-trinity-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Run Trinity to assemble a transcriptome using a directory of FASTQ files."
  echo
  echo "## Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING     Input directory with FASTQ files"
  echo "## -o STRING     Output directory"
  echo "                 NOTE: The output directory needs to include 'trinity' in its name"
  echo
  echo "## Other options:"
  echo "## -h            Print this help message and exit"
  echo
  echo "## Example: $0 -i data/fastq/ -o results/trinity"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
indir=""
outdir=""

## Get parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script trinity.sh"
date
echo

## Check parameter values
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-i) $indir does not exist" >&2 && exit 1


# SOFTWARE ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
CONDA_ENV=/users/PAS0471/jelmer/miniconda3/envs/trinity-env
COLLECTL_SCRIPT="$CONDA_ENV"/opt/trinity-2.12.0/trinity-plugins/COLLECTL/examine_resource_usage_profiling.pl
source activate "$CONDA_ENV"


# OTHER SETUP ------------------------------------------------------------------ 
## Bash strict mode
set -euo pipefail

## Comma-delimited list of FASTQ files:
R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')

## Other parameters
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G

## Report
echo "## Input dir:                            $indir"
echo "## Output dir:                           $outdir"
echo
echo "## List of R1 (forward) FASTQ files:     $R1_list"
echo "## List of R2 (reverse) FASTQ files:     $R2_list"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN TRINITY ---------------------------------------------------------------
echo "## Starting Trinity run..."
Trinity --seqType fq \
        --left "$R1_list" \
        --right "$R2_list" \
        --SS_lib_type RF \
        --output "$outdir" \
        --max_memory "$mem_gb" \
        --CPU "$SLURM_CPUS_ON_NODE" \
        --verbose

#--monitoring \

## Check resource usage - https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Runtime-Profiling
#mv collectl "$outdir"
#"$COLLECTL_SCRIPT" "$outdir"/collectl


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script trinity.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
