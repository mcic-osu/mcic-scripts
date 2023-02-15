#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-qcat_demux-%j.out

# SETUP ------------------------------------------------------------------------
## Load software 
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/qcat-env

## Bash strict settings
set -ueo pipefail

## Command-line args
infile=$1
outdir=$2

## Input checks
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments; you provided $#" && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input infile $infile does not exist" && exit 1

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"

## Report
echo "## Starting script qcat_demux.sh..."
date
echo
echo "## Input FASTQ file:             $infile"
echo "## Output dir:                   $outdir"
echo -e "--------------------\n"


# RUN --------------------------------------------------------------------------
qcat -f "$infile" -b "$outdir"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing files in output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script qcat_demux.sh"
date


#? https://github.com/nanoporetech/qcat