#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --job-name=prokka
#SBATCH --output=slurm-prokka-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Prokka to annotate a prokaryotic genome assembly"
  echo
  echo "Syntax: $0 -i <input-file> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Input file, a nucleotide FASTA file"
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to Prokka"
  echo "    -g STRING         Genus name"
  echo "    -s STRING         Species name"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:     $0 -i results/spades/assembly.fasta -o results/prokka"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Prokka documentation: https://github.com/tseemann/prokka"
  echo "Prokka paper: https://pubmed.ncbi.nlm.nih.gov/24642063/"
  echo
}

## Option defaults
infile=""
outdir=""
genus=""
species=""
more_args=""

## Parse command-line options
while getopts ':i:o:g:s:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    g) genus="$OPTARG" ;;
    s) species="$OPTARG" ;;
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


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/prokka-1.14.6

## Bash strict settings
set -euo pipefail

## Define output files etc
sampleID=$(basename "$infile" .fasta)

## Build genus and species arguments
if [[ "$genus" != "" ]]; then
    genus_arg="--genus $genus"
else
    genus_arg=""
fi

if [[ "$species" != "" ]]; then
    species_arg="--species $species"
else
    species_arg=""
fi

## Report
echo "## Starting script prokka.sh"
date
echo
echo "## Input file:                           $infile"
echo "## Output dir:                           $outdir"
[[ "$genus" != "" ]] && echo "## Genus:                                $genus"
[[ "$species" != "" ]] && echo "## Species:                              $species"
[[ $more_args != "" ]] && echo "## Other arguments to pass to prokka:    $more_args"
echo -e "--------------------\n"

## Make output dir
mkdir -p "$outdir"


# RUN prokka --------------------------------------------------------------------
echo -e "\n## Now running prokka...."
prokka \
    --outdir "$outdir" \
    --prefix "$sampleID" \
    --cpus "$SLURM_CPUS_PER_TASK"  \
    --force $genus_arg $species_arg $more_args \
    "$infile"

# --strain
# --usegenus  ?
# --addgenes

## Remove DNA sequences from GFF file
echo -e "\n## Now removing DNA sequences from GFF file..."
mv "$outdir"/"$sampleID".gff "$outdir"/"$sampleID"_withseqs.gff
sed '/^##FASTA/Q' "$outdir"/"$sampleID"_withseqs.gff > "$outdir"/"$sampleID".gff


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"/"$sampleID"*
echo -e "\n## Done with script prokka.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
