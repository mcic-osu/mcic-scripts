#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=evigene
#SBATCH --output=slurm-evigene-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run EviGene to filter a transcriptome assembly."
  echo
  echo "Syntax: $0 -i <transcriptome-FASTA> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Transcriptome assembly FASTA file"
  echo "                      NOTE: Concatenate multiple assemblies FASTAs into a single one prior to running this script."
  echo "    -o STRING         Output dir"
  echo "                      NOTE: Any files present in this dir will be removed prior to running EviGene!"
  echo
  echo "Other options:"
  echo "    -m INTEGER        Minimum CDS size            [default: 350]"
  echo "    -a STRING         Other argument(s) to pass to EviGene"
  echo "    -h                Print this help message"
  echo
  echo "Example:              $0 -i results/merged_assembly.fasta -o results/evigene"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "EviGene documentation: http://arthropods.eugenes.org/EvidentialGene/evigene/"
  echo "EviGene output:        http://arthropods.eugenes.org/EvidentialGene/evigene/docs/EvigeneR/evigene4_outputs_brief.txt"
  echo "EviGene update note:   http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_update2020march.html"
  echo "EviGene paper:         https://www.biorxiv.org/content/10.1101/829184v1"
  echo "Evigene repo:          https://sourceforge.net/projects/evidentialgene"
  echo
}

## Option defaults
infile=""
outdir=""
min_cds=350
more_args=""

## Parse command-line options
while getopts ':i:o:m:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    m) min_cds="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Check input
[[ "$infile" = "" ]] && echo "## ERROR: Please specify an input dir with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file (-i) $infile does not exist" >&2 && exit 1


# SET-UP -----------------------------------------------------------------------
## Software
CONDA_ENV_DIR=/fs/project/PAS0471/jelmer/conda/evigene
module load python/3.6-conda5.2
source activate $CONDA_ENV_DIR

## Bash strict mode
set -euo pipefail

## Report
echo
echo "## Starting script evigene.sh"
date
echo
echo "## Input FASTA file:                     $infile"
echo "## Output dir:                           $outdir"
echo "## Minimum CDS size:                     $min_cds"
[[ "$more_args" != "" ]] && echo "## Other arguments to pass to evigene:   $more_args"
echo -e "--------------------\n"

## Remove output dir if it already exists (EviGene behaves weirdly when there are files present!)
[[ "$outdir" != "" && -d "$outdir" ]] && rm -r "${outdir:?}"

## Make output dir
mkdir -p "$outdir"

## Copy input file to outdir
echo "## Copying input FASTA to output dir..."
infile_base=$(basename "$infile")
cp -v "$infile" "$outdir"


# RUN EVIGENE ------------------------------------------------------------------
cd "$outdir" || exit

echo -e "\n## Now running Evigene..."
$CONDA_ENV_DIR/bin/prot/tr2aacds.pl \
    -debug \
    -MINCDS "$min_cds" \
    -NCPU "$SLURM_CPUS_PER_TASK" \
    -MAXMEM "$SLURM_MEM_PER_NODE" \
    -log \
    -tidyup $more_args \
    -mrnaseq "$infile_base"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh
echo -e "\n## Done with script evigene.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
