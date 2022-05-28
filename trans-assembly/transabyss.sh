#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=transabyss
#SBATCH --output=slurm-transabyss-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Trans-ABySS to create a transcriptome assembly."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Genome (nucleotide) FASTA file"
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to Trans-ABySS"
  echo "    -h                Print this help message"
  echo
  echo "Example: $0 -i data/fastq -o results/transabyss"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Trans-ABySS documentation: https://github.com/bcgsc/transabyss/blob/master/TUTORIAL.md"
  echo
}

## Option defaults
indir=""
outdir=""
assembly_ID=""
kmer_size=32
minlen=350
more_args=""

## Parse command-line options
while getopts ':i:I:o:k:l:a:h' flag; do
  case "${flag}" in
    i) indir="$OPTARG" ;;
    I) assembly_ID="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    k) kmer_size="$OPTARG" ;;
    l) minlen="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Check input
[[ "$indir" = "" ]] && echo "## ERROR: Please specify an input dir with -i" >&2 && exit 1
[[ "$assembly_ID" = "" ]] && echo "## ERROR: Please specify an assembly ID with -I" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-i) $indir does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/transabyss-2.0.1

## Bash strict mode
set -euo pipefail

## Make output dir
mkdir -p "$outdir"

## Report
echo
echo "## Starting script transabyss.sh"
date
echo
echo "## Input dir:                            $indir"
echo "## Output dir:                           $outdir"
echo "## Assembly ID:                          $assembly_ID"
echo
echo "## Kmer size:                            $kmer_size"
echo "## Min contig length:                    $minlen"
[[ "$more_args" != "" ]] &&echo "## Other arguments for Trans-ABySS:       $more_args"
echo -e "--------------------\n"


# RUN SOAPDENOVO ---------------------------------------------------------------
echo -e "\n## Now running Trans-ABySS..."
transabyss \
    --pe "$indir"/* \
    --SS \
    --kmer "$kmer_size" \
    --length "$minlen" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --outdir "$outdir" \
    --name "$assembly_ID" $more_args

#--se trimmed_unpaired*/* \


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script transabyss.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
