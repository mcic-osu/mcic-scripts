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
  echo "$0: Normalize RNAseq reads prior to assembly with ORNA"
  echo
  echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "-i STRING        Input directory with FASTQ files"
  echo "-o STRING        Output directory"
  echo
  echo "Other options:"
  echo "-h               Print this help message and exit"
  echo
  echo "Example: $0 -i data/fastq/ -o results/orna"
  echo "To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "ORNA docs:       https://github.com/SchulzLab/ORNA"
  echo "ORNA paper:      https://www.nature.com/articles/s41598-019-41502-9"
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

## Report
echo -e "\n## Starting script orna.sh"
date
echo

## Check parameter values
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file R1 (-i) $R1_in does not exist" >&2 && exit 1


# SOFTWARE ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/orna-2.0


# OTHER SETUP ------------------------------------------------------------------ 
## Bash strict mode
set -euo pipefail

## Determine name of R2 file
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}

## Determine output prefix
R1_basename=$(basename "$R1_in" .fastq.gz)
sampleID=${R1_basename/"$R1_suffix"/}
prefix="$outdir"/"$sampleID"

## Report
echo "## R1 input file:              $R1_in"
echo "## R2 input file:              $R2_in"
echo "## Output prefix:              $prefix"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN KHMER ---------------------------------------------------------------
echo "## Starting normalization ..."
ORNA \
    -pair1 "$R1_in" -pair2 "$R2_in" \
    -output "$prefix" \
    -type fastq \
    -nb-cores "$SLURM_CPUS_PER_TASK"

## Move output files
R2_basename=$(basename "$R2_in" .fastq.gz)
mv ./*"$R2_basename"*h5 "$outdir"
gzip -c "$prefix"_1.fq > "$prefix"_R1.fastq.gz
gzip -c "$prefix"_2.fq > "$prefix"_R2.fastq.gz


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script orna.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
