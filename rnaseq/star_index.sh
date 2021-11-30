#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=STAR_index
#SBATCH --output=slurm-STAR-index-%j.out


# SETUP ---------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env

## Strict bash settings
set -euo pipefail

## Help function
Help() {
  echo
  echo "## $0: Index a reference genome FASTA file with STAR."
  echo
  echo "## Syntax: $0 -i <R1 FASTQ input file> -o <output dir> [ -s <index size> ] [-sh]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   R1 FASTQ input file (REQUIRED)"
  echo "## -o STR   Output dir for index files (REQUIRED)"
  echo "## -s INT   Index size (default: 13)"
  echo "## Example: $0 -i refdata/my_genome.fa -o refdata/star_index -s 10"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
ref_fa=""
index_dir=""
index_size=13

## Parse command-line options
while getopts ':i:o:s:h' flag; do
  case "${flag}" in
  i) ref_fa="$OPTARG" ;;
  o) index_dir="$OPTARG" ;;
  s) index_size="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Report
echo "## Starting script star_index.sh"
date 
echo "## Genome FASTA file (input): $ref_fa"
echo "## Genome index dir (output): $index_dir"
echo "## Index size (genomeSAindexNbases): $index_size"
echo -e "-----------------------\n"

## Make output dir if needed
mkdir -p "$index_dir"

## STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $ref_fa = *gz ]]; then
    echo "## Unzipping gzipped FASTA file"
    gunzip "$ref_fa"
    ref_fa=${ref_fa/.gz/}
    echo "## Genome FASTA file: $ref_fa"
fi


# RUN STAR ---------------------------------------------------------------------
echo "## Running STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
     --runMode genomeGenerate \
     --genomeDir "$index_dir" \
     --genomeFastaFiles "$ref_fa" \
     --genomeSAindexNbases "$index_size" \
     --sjdbOverhang ReadLength-1


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$index_dir"
echo -e "\n## Done with script star_index.sh"
date
