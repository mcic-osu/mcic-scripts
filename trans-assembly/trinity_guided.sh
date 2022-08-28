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
  echo "$0: Run Trinity to assemble a genome-guided transcriptome using a BAM file."
  echo
  echo "Syntax: $0 -i <input-BAM> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "   -i STRING       Input BAM file"
  echo "   -o STRING       Output directory"
  echo "                   NOTE: The output directory needs to include 'trinity' in its name"
  echo
  echo "Other options:"
  echo "   -I INTEGER      Max. intron size                     [default: 10000]"
  echo "   -h              Print this help message and exit"
  echo
  echo "Example: $0 -i results/star/my.bam -o results/trinity"
  echo "To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
bam=""
outdir=""
max_intron_size=10000

## Get parameter values
while getopts ':i:I:o:h' flag; do
    case "${flag}" in
    i) bam="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    I) max_intron_size="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script trinity_guided.sh"
date
echo

## Check parameter values
[[ ! -f "$bam" ]] && echo "## ERROR: Input BAM file (-i) $bam does not exist" >&2 && exit 1


# SOFTWARE ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
CONDA_ENV=/users/PAS0471/jelmer/miniconda3/envs/trinity-env
COLLECTL_SCRIPT="$CONDA_ENV"/opt/trinity-2.12.0/trinity-plugins/COLLECTL/examine_resource_usage_profiling.pl
source activate "$CONDA_ENV"


# OTHER SETUP ------------------------------------------------------------------ 
## Bash strict mode
set -euo pipefail

## Other parameters
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G

## Report
echo "## Input BAM file:                       $bam"
echo "## Output dir:                           $outdir"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN TRINITY ---------------------------------------------------------------
echo "## Starting Trinity run..."
Trinity --genome_guided_bam "$bam" \
        --genome_guided_max_intron "$max_intron_size" \
        --output "$outdir" \
        --max_memory "$mem_gb" \
        --CPU "$SLURM_CPUS_ON_NODE" \
        --monitoring \
        --verbose

## Check resource usage - https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Runtime-Profiling
mv collectl "$outdir"
"$COLLECTL_SCRIPT" "$outdir"/collectl


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script trinity_guided.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo

