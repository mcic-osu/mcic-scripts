#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --job-name=braker
#SBATCH --output=slurm-braker-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Braker2 to annotate a genome."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -o <output-dir> -d <OrthoDB-FASTA> ..."
  echo
  echo "Required options:"
  echo "-i STRING         Genome (nucleotide) FASTA file"
  echo "-o STRING         Output dir"
  echo "-d STRING         OrthoDB (protein) FASTA file"
  echo "                  For info on how to create this file: https://github.com/gatech-genemark/ProtHint#protein-database-preparation"
  echo
  echo "Other options:"
  echo "-a STRING         Other argument(s) to pass to Braker2"
  echo "-h                Print this help message"
  echo
  echo "Example: $0 -i my_genome.fa -o results/braker -d odb_prots.fa"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Braker2 documentation: https://github.com/Gaius-Augustus/BRAKER"
  echo "Braker2 paper: https://academic.oup.com/nargab/article/3/1/lqaa108/6066535"
  echo
}

## Option defaults
genome_fa=""
orthodb_protein_fa=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:d:ah' flag; do
  case "${flag}" in
    i) genome_fa="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) orthodb_protein_fa="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Report
echo "## Starting script braker.sh"
date
echo

## Check input
[[ ! -f "$genome_fa" ]] && echo "## ERROR: Input file (-i) $genome_fa does not exist" >&2 && exit 1
[[ ! -f "$orthodb_protein_fa" ]] && echo "## ERROR: OrthoDB protein file (-d) $orthodb_protein_fa does not exist" >&2 && exit 1


# LOAD SOFTWARE ----------------------------------------------------------------
## Braker2 conda env which contains everything except GeneMark-EX and ProtHint 
conda activate /fs/project/PAS0471/jelmer/conda/braker2-env

## GeneMark-EX
# See https://github.com/Gaius-Augustus/BRAKER#genemark-ex
GENEMARK_BASEDIR=/fs/project/PAS0471/jelmer/software/genemark-ex
export GENEMARK_PATH="$GENEMARK_BASEDIR"/gmes_linux_64_4
cp "$GENEMARK_BASEDIR"/gm_key_64 ~/.gm_key

# See https://github.com/Gaius-Augustus/BRAKER#prothint
export PROTHINT_PATH=/fs/project/PAS0471/jelmer/software/ProtHint


# OTHER SETUP ------------------------------------------------------------------
## Bash strict mode
set -euo pipefail

## Report
echo
echo "## Genome FASTA file:                    $genome_fa"
echo "## OrthoDB protein FASTA:                $orthodb_protein_fa"
echo "## Output dir:                           $outdir"
echo "## Other arguments to pass to Braker:    $more_args"
echo -e "--------------------\n"

## If needed, make dirs absolute because we have to move into the outdir
[[ $genome_fa =~ ^/ ]] && genome_fa="$PWD"/"$genome_fa"
[[ $orthodb_protein_fa =~ ^/ ]] && orthodb_protein_fa="$PWD"/"$orthodb_protein_fa"


# RUN BRAKER2 -----------------------------------------------------------------
echo "## Now running Braker2..."

cd "$outdir" || exit 1

braker.pl \
    --genome="$genome_fa" \
    --prot_seq="$orthodb_protein_fa" \
    --softmasking \
    $more_args \
    --cores="$SLURM_CPUS_PER_TASK"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script braker.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
