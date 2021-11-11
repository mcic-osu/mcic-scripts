#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=STAR_align
#SBATCH --output=slurm-STAR-align-%j.out


# SETUP ---------------------------------------------------------------------
## Load software
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source ~/.bashrc
source activate star-env
conda activate --stack samtools-env

## Bash strict mode
set -euo pipefail

## Help
Help() {
  echo
  echo "## $0: Index a reference genome FASTA file with STAR."
  echo
  echo "## Syntax: $0 -i <R1 FASTQ input file> -o <BAM output dir> -r <ref genome index dir> [ -m <max multi-map> ] [-h]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   R1 FASTQ input file (REQUIRED; note that the name of the R2 file will be inferred by the script.)"
  echo "## -o STR   BAM output dir (REQUIRED)"
  echo "## -r STR   Reference index dir (REQUIRED)"
  echo "## -m INT   Max. number of locations a read can map to, before being considered unmapped (default: 10)"
  echo "## Example: $0 -i data/fastq/S01_L001_R1.fastq.gz -o results/star -r refdata/star_index"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
R1_in=""
bam_dir=""
index_dir=""
max_map=10

# Get command-line options:
while getopts ':i:o:r:m:h' flag; do
  case "${flag}" in
  i) R1_in="$OPTARG" ;;
  o) bam_dir="$OPTARG" ;;
  r) index_dir="$OPTARG" ;;
  m) max_map="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Process args
fq_R2=${R1_in/_R1_/_R2_}
sampleID=$(basename "$R1_in" | sed 's/_R1_.*//')
dir_flagstat="$bam_dir"/flagstat
flagstat_out="$dir_flagstat"/"$sampleID"_flagstat.txt

## Report
echo "## Starting script star_align.sh"
date
echo "## R1 FASTQ file (input): $R1_in"
echo "## R2 FASTQ file (input - inferred): $fq_R2"
echo "## Sample ID (inferred): $sampleID"
echo "## Genome index dir (input): $index_dir"
echo "## BAM dir (output): $bam_dir"
echo "## Max nr of alignments for a read (setting): $max_map"  # If exceeded, read is considered unmapped
echo -e "------------------------\n"

## Check inputs
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file $R1_in does not exist" && exit 1
[[ ! -f "$fq_R2" ]] && echo "## ERROR: Input file $fq_R2 does not exist" && exit 1
[[ ! -d "$index_dir" ]] && echo "## ERROR: Input dir $index_dir does not exist" && exit 1

## Create output dirs if needed
mkdir -p "$bam_dir"
mkdir -p "$dir_flagstat"

# ALIGN ------------------------------------------------------------------------
echo "## Running STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
   --genomeDir "$index_dir" \
   --readFilesIn "$R1_in" "$fq_R2" \
   --readFilesCommand zcat \
   --outFileNamePrefix "$bam_dir/$sampleID"_ \
   --outFilterMultimapNmax $max_map \
   --outSAMtype BAM SortedByCoordinate \
   --outReadsUnmapped Fastx \
   --alignIntronMin 5 --alignIntronMax 350000

## Move logfiles
mkdir -p "$bam_dir"/star_logs  
mv "$bam_dir"/"$sampleID"*out "$bam_dir"/"$sampleID"*tab "$bam_dir"/star_logs/


# RUN FLAGSTAT -----------------------------------------------------------------
echo "## Running samtools flagstat...."
samtools flagstat "$bam_dir/$sampleID"_*bam > "$flagstat_out"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n\n----------------------"
echo "## Listing output BAM file:"
ls -lh "$bam_dir/$sampleID"_*bam
echo "## Listing flagstat file:"
ls -lh "$flagstat_out"

echo -e "\n## Done with script star_align.sh"
date
