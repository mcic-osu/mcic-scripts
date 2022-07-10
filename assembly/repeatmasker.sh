#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --job-name=repeatmasker
#SBATCH --output=slurm-repeatmasker-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run RepeatMasker to repeat-mask a genome."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -o <output-dir> -s <species> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Genome (nucleotide) FASTA file"
  echo "    -o STRING         Output dir"
  echo "    -s STRING         Species or taxonomic group name"
  echo "                      To check which species/groups are available, run, e.g:"
  echo "                      /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'oomycetes'"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to RepeatMasker"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i my_genome.fa -o results/repeatmasker -s stramenopiles"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "RepeatMasker documentation: https://www.repeatmasker.org/"
  echo
}

## Option defaults
genome_fa=""
species=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:s:a:h' flag; do
  case "${flag}" in
    i) genome_fa="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    s) species="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

# SETUP ------------------------------------------------------------------------
## Check input
[[ ! -f "$genome_fa" ]] && echo "## ERROR: Input file (-i) $genome_fa does not exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please provide an output dir with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1

## Bash script settings
set -euo pipefail

## Species arg
if [[ "$species" != "" ]]; then
    species_arg="-species $species"
else
    species_arg=""
fi

## Report
echo
echo "## Starting script repeatmasker.sh"
date
echo
echo "## Input file (genome FASTA):                  $genome_fa"
echo "## Output dir:                                 $outdir"
[[ $species_arg != "" ]] && echo "## Species name:                               $species"
[[ $more_args != "" ]] && echo "## Other arguments to pass to RepeatMasker:    $more_args"
echo -e "--------------------\n"

## Make output dir
mkdir -p "$outdir"

# RUN REPEATMASKER -------------------------------------------------------------
echo "## Now runnning RepeatMasker..."
RepeatMasker \
    -dir "$outdir" \
    $species_arg $more_args "$genome_fa"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script repeatmasker.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo

## To check available species, e.g:
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'oomycetes'
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'stramenopiles'
# /fs/project/PAS0471/jelmer/conda/repeatmasker-4.1.2.p1/share/RepeatMasker/famdb.py names 'phytophthora'
#>67593 Phytophthora megasperma f. sp. glycinea (includes), Phytophthora sojae (scientific name)