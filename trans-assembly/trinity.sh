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

## Report
echo -e "\n## Starting script trinity.sh"
date

## Command-line args
indir="$1"
outdir="$2"      # NOTE: Output dir needs to include "trinity" in the name according to the Trinity docs

## Check input
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir $indir does not exist" && exit 1

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/trinity-env
COLLECTL_SCRIPT=/users/PAS0471/jelmer/miniconda3/envs/trinity-env/opt/trinity-2.12.0/trinity-plugins/COLLECTL/examine_resource_usage_profiling.pl

## Bash strict mode
set -euo pipefail

## Comma-delimited list of FASTQ files:
R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')

## Other parameters
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G

## Report
echo "## All args: $*"
echo
echo "## Input dir:         $indir"
echo "## Output dir:        $outdir"
echo
echo "## List of R1 files:  $R1_list"
echo "## List of R2 files:  $R2_list"
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
        --monitoring \
        --verbose

## Check resource usage - https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Runtime-Profiling
mv collectl "$outdir"
"$COLLECTL_SCRIPT" "$outdir"/collectl


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script trinity.sh"
date
echo
