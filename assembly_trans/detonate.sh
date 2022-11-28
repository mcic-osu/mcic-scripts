#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=detonate
#SBATCH --output=slurm-detonate-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Detonate to evaluate a transcriptome assembly."
  echo
  echo "Syntax: $0 -i <input-FASTA> -I <input-FASTQ-dir> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Input FASTA file (transcriptome assembly)"
  echo "    -I STRING         Input dir with FASTQ files"
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -l INTEGER        FASTQ read length          [default: 150]"
  echo "    -a STRING         Other argument(s) to pass to Detonate"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i results/assembly/contigs.fasta -I data/fastq -o results/detonate"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Detonate website: http://deweylab.biostat.wisc.edu/detonate/vignette.html"
  echo "Detonate paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4298084/"
  echo
}

## Option defaults
fa_in=""
fq_indir=""
outdir=""
more_args=""
readlen=150

## Parse command-line options
while getopts ':i:I:o:l:a:h' flag; do
  case "${flag}" in
    i) fa_in="$OPTARG" ;;
    I) fq_indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    l) readlen="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Other parameters
R1_LIST=$(echo "$fq_indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_LIST=$(echo "$fq_indir"/*R2*fastq.gz | sed 's/ /,/g')
N_CORES=$SLURM_CPUS_PER_TASK
TRANS_ID=$(basename "$fa_in")
PREFIX="$outdir"/"${TRANS_ID%.fa*}"

## Check input
[[ "$fa_in" = "" ]] && echo "## ERROR: Please specify an input FASTQ file with -i" >&2 && exit 1
[[ "$fq_indir" = "" ]] && echo "## ERROR: Please specify a FASTQ input dir with -I" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file (-i) $fa_in does not exist" >&2 && exit 1
[[ ! -d "$fq_indir" ]] && echo "## ERROR: Input dir (-I) $fq_indir does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/detonate-env

## Bash strict settings
set -euo pipefail

## Report
echo
echo "## Starting script detonate.sh"
date
echo "## Input FASTA file:              $fa_in"
echo "## Input dir with FASTQ files:    $fq_indir"
echo "## Output dir:                    $outdir"
echo "## Read length:                   $readlen"
echo "## Full output prefix:            $PREFIX"
echo "## Number of cores:               $N_CORES"
echo
echo "## List of R1 files:              $R1_LIST"
echo
echo "## List of R2 files:              $R2_LIST"
echo -e "-------------------------------\n"


# RUN DETONATE -----------------------------------------------------------------
## Make output dir
mkdir -p "$outdir"/logs

echo -e "\n## Running rsem-eval-calculate-score..."
rsem-eval-calculate-score \
    --paired-end \
    "$R1_LIST" "$R2_LIST" \
    "$fa_in" \
    "$PREFIX" \
    "$readlen" \
    -p "$N_CORES" \
    $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script detonate.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo


# INFO -------------------------------------------------------------------------
#> RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic
#> model that depends only on an assembly and the RNA-Seq reads used for its construction.
#> Unlike N50, RSEM-EVAL combines multiple factors, including the compactness of an assembly
#> and the support of the assembly from the RNA-Seq data, into a single, statistically-principled evaluation score.
#> This score can be used to select a best assembler, optimize an assembler's parameters, and guide new assembler design as an objective function.
#> In addition, for each contig within an assembly, RSEM-EVAL provides a score that
#> assesses how well that contig is supported by the RNA-Seq data and can be used to filter unnecessary contigs.
