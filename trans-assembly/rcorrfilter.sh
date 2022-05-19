#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=rcorrfilter
#SBATCH --output=slurm-rcorrfilter-%j.out

# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Filter an rcorrector-processed pair of FASTQ files."
  echo
  echo "## Syntax: $0 -i <input-FASTA> -o <output-dir> -a <gff-file> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING     R1 input FASTQ file as output by rcorrector.sh (name of R2 file will be inferred)"
  echo "## -o STRING     Output directory for corrected FASTQ FILES"
  echo
  echo "## Other options:"
  echo "## -h            Print this help message and exit"
  echo
  echo "## Example: $0 -i results/rcorr/A1_R1.fastq.gz -o results/readfilt"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameters
R1=""
outdir=""

## Get parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) R1="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script rcorrfilter.sh"
date
echo

## Test parameter values
[[ ! -f "$R1" ]] && echo "## ERROR: Input file (-i) $R1 does not exist" >&2 && exit 1


# SOFTWARE ---------------------------------------------------------------------
## Software and scripts
module load python/3.6-conda5.2
CONDA_ENV=/users/PAS0471/jelmer/miniconda3/envs/rcorrector-env
FILTER_SCRIPT=$CONDA_ENV/bin/FilterUncorrectabledPEfastq.py
source activate "$CONDA_ENV"
#? The script `FilterUncorrectabledPEfastq.py` from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
#? has been added to the rcorrector Conda environment


# OTHER SETUP ------------------------------------------------------------------
## Bash strict mode
set -euo pipefail

## Process parameter values
R2=${R1/_R1/_R2}
sample_id=$(echo "$(basename $R1)" | sed 's/_R1.*//')

R1_out="$outdir"/$(basename "$R1" .cor.fq.gz).fastq
R2_out="$outdir"/$(basename "$R2" .cor.fq.gz).fastq

## Report
echo "## R1 input file:               $R1"
echo "## R2 input file:               $R2"
echo "## Output dir:                  $outdir"
echo "## R1 output file:              $R1_out"
echo "## R2 output file:              $R2_out"
echo "## Sample ID:                   $sample_id"
echo -e "---------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN PYTHON SCRIPT ------------------------------------------------------------
echo "## Running read-filter script..."
python2 "$FILTER_SCRIPT" -s "$sample_id" -1 "$R1" -2 "$R2"


# WRAP-UP ----------------------------------------------------------------------
## Move output files
echo -e "\n## Moving and gzipping output files..."
echo unfixrm*"$sample_id"*R1*cor.fq
mv -v unfixrm_"$sample_id"*R1*cor.fq "$R1_out"
gzip -f "$R1_out"

echo unfixrm*"$sample_id"*R2*cor.fq
mv -v unfixrm_"$sample_id"*R2*cor.fq "$R2_out"
gzip -f "$R2_out"

mv rmunfixable_"$sample_id".log "$outdir"

## Report
echo -e "\n-------------------------------"
echo "## Listing output files:"
ls -lh "$R1_out".gz "$R2_out".gz

echo -e "\n## Done with script rcorrfilter.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
