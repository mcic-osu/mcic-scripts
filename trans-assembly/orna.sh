#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=orna
#SBATCH --output=slurm-orna-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run ORNA to normalize paired-end RNAseq reads prior to assembly"
  echo
  echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i FILE             Input R1 FASTQ file (R2 filename will be inferred)"
  echo "    -o DIR              Output directory"
  echo
  echo "Other options:"
  echo "    -h                  Print this help message and exit"
  echo
  echo "Example:                $0 -i data/fastq/ -o results/orna"
  echo "To submit this script to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "ORNA docs:  https://github.com/SchulzLab/ORNA"
  echo "ORNA paper: https://www.nature.com/articles/s41598-019-41502-9"
  echo
}

## Default parameter values
R1_in=""
outdir=""

## Get parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) R1_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Check parameter values
[[ "$R1_in" = "" ]] && echo "## ERROR: Please specify an input file with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output directory with -o" >&2 && exit 1
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file R1 (-i) $R1_in does not exist" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/orna-2.0

## Bash strict mode
set -euo pipefail

## Determine name of R2 file
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
[[ "$R1_in" = "$R2_in" ]] && echo "## ERROR: Input file R1 is the same as R2" >&2 && exit 1

## Determine output prefix
R1_basename=$(basename "$R1_in" | sed -E 's/.fa?s?t?q.gz//')
sampleID=${R1_basename/"$R1_suffix"/}

## If needed, make paths absolute because we have to move into the outdir
[[ ! $R1_in =~ ^/ ]] && R1_in="$PWD"/"$R1_in"
[[ ! $R2_in =~ ^/ ]] && R2_in="$PWD"/"$R2_in"

## Report
echo
echo "## Starting script orna.sh"
date
echo
echo "## R1 input file:              $R1_in"
echo "## R2 input file:              $R2_in"
echo "## Sample ID:                  $sampleID"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN ORNA ---------------------------------------------------------------------
cd "$outdir" || exit

echo "## Starting normalization ..."
ORNA \
    -pair1 "$R1_in" \
    -pair2 "$R2_in" \
    -output "$sampleID" \
    -type fastq \
    -nb-cores "$SLURM_CPUS_PER_TASK"

## Rename and gzip output files
echo "## Compressing output FASTQ files..."
gzip -cv "$sampleID"_1.fq > "$sampleID"_R1.fastq.gz && rm "$sampleID"_1.fq
gzip -cv "$sampleID"_2.fq > "$sampleID"_R2.fastq.gz && rm "$sampleID"_2.fq
rm "$sampleID"*h5


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$sampleID"*
echo -e "\n## Done with script orna.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
