#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm-trimgalore-%j.out


# SETUP ------------------------------------------------------------------------
## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/trimgalore-env

## Bash strict settings
set -euo pipefail

## Help
Help() {
  echo
  echo "## $0: Run TrimGalore for a FASTQ file."
  echo
  echo "## Syntax: $0 -i <R1 FASTQ input file> -o <FASTQ output dir> -O <FastQC output dir> -q <min seq qual>  -q <min seq len> [-sh]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   R1 FASTQ input file (REQUIRED)"
  echo "## -o STR   FASTQ output dir (REQUIRED)"
  echo "## -O STR   FastQC results output dir (REQUIRED)"
  echo "## -q INT   Quality trimming threshold (default: 20)"
  echo "## -l INT   Minimum read length (default: 20)"
  echo "## -s       Input is single-end (default: paired-end)"
  echo "## Example: $0 -i data/fastq/S01_L001_R1.fastq.gz -o results/trimgalore -O results/fastqc"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
qual=0                  # => no quality trimming
len=20                  # => 20 is also the TrimGalore default
single_end=fastq_file   # => paired-end by defaults

# Get command-line options:
while getopts ':i:o:O:q:l:sh' flag; do
  case "${flag}" in
  i) R1_in="$OPTARG" ;;
  o) outdir_trim="$OPTARG" ;;
  O) outdir_fastqc="$OPTARG" ;;
  q) qual="$OPTARG" ;;
  l) len="$OPTARG" ;;
  s) single_end=true ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Process variables and args
if [ "$single_end" != "true" ]; then
    R2_in=${R1_in/_R1_/_R2_}
    R2_id=$(basename "$R2_in" .fastq.gz)
    input_arg="--paired $R1_in $R2_in"
else
    input_arg="$R1_in"
fi

R1_id=$(basename "$R1_in" .fastq.gz)

n_threads="$SLURM_CPUS_PER_TASK"

## Make output dirs
mkdir -p "$outdir_trim"
mkdir -p "$outdir_fastqc"

## Report
echo "## Starting script trimgalore.sh"
date
echo "## R1 input file:                 $R1_in"
echo "## Output dir - trimmed FASTQs:   $outdir_trim"
echo "## Output dir - FastQC:           $outdir_fastqc"
echo "## Sequence quality threshold:    $qual"
echo "## Minimum sequence length:       $len"
[[ "$single_end" != "true" ]] && echo "## R2 input file:                 $R2_in"
echo -e "---------------------------\n\n"


# RUN TRIMGALORE AND PROCESS OUTPUT --------------------------------------------
## Run Trim-Galore
trim_galore \
    --quality "$qual" --length "$len" \
    --gzip -j "$n_threads" \
    --output_dir "$outdir_trim" \
    --fastqc --fastqc_args "-t $n_threads --outdir $outdir_fastqc" \
    $input_arg

## Move FASTQ files
if [ "$single_end" != "true" ]; then
    mv "$outdir_trim"/"$R1_id"_val_1.fq.gz "$outdir_trim"/"$R1_id".fastq.gz
    mv "$outdir_trim"/"$R2_id"_val_2.fq.gz "$outdir_trim"/"$R2_id".fastq.gz
else
    mv "$outdir_trim"/"$R1_id"_trimmed.fq.gz "$outdir_trim"/"$R1_id".fastq.gz
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir_trim"/"$R1_id".fastq.gz "$outdir_trim"/"$R2_id".fastq.gz

echo -e "\n## Done with script trimgalore.sh"
date