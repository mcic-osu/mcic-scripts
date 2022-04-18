#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --output=slurm-extract_kraken_reads-%j.out

# ARGS AND PARAMETERS ----------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Extract reads assigned to a given taxon ID by Kraken2"
    echo
    echo "Syntax: $0 -i <input-R1-FASTQ> -I <kraken-output-file> -t <taxon-ID(s) -o <output dir> ..."
    echo
    echo "Required options:"
    echo "  -i STRING     Input R1 (forward reads) sequence file (name of R2 will be inferred)"
    echo "  -I STRING     Kraken 'main' output file with per-read taxon. assignments"
    echo "                The FASTQ and Kraken output file should be for the same sample"
    echo "  -t STRING     NCBI Taxon IDs -- if giving multiple, separate by spaces and quote, e.g. '699189 1111709'"
    echo "  -o STRING     Output directory (will be created if needed)"
    echo
    echo "Other options:"
    echo "  -h            Print this help message and exit"
    echo
    echo "Example command:"
    echo "$0 -i data/fastq/A1_R1.fastq.gz -I results/kraken/A1_main.txt -t '699189 1111709' -o results/kraken/reads"
    echo
    echo "To submit to the OSC queue, preface with sbatch:"
    echo "sbatch $0 ..."
    echo
    echo "SLURM parameters in script: '--account=PAS0471 --time=1:00:00 --output=slurm-extract_kraken_reads-%j.out"
    echo
    echo "Default SLURM parameters can be overridden when submitting the script, e.g.:"
    echo "sbatch -t 15 $0 ...      (override default time reservation of 3 hours, use 15 minutes instead)"
    echo
    echo "Krakentools documentation: https://github.com/jenniferlu717/KrakenTools"
    echo
}

## Option defaults
R1_in=""
outdir=""
kraken_main=""
taxids=""

## Parse options
while getopts ':i:o:I:t:h' flag; do
    case "${flag}" in
    i) R1_in="$OPTARG" ;;
    I) kraken_main="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    t) taxids="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Report
echo "## Starting script extract_kraken_reads.sh"
date
echo

## Check input
[[ "$R1_in"  = "" ]] && echo "## ERROR: Please provide R1 input FASTQ file with -i flag" && exit 1
[[ "$outdir"  = "" ]] && echo "## ERROR: Please provide output dir with -o flag" && exit 1
[[ "$taxids"  = "" ]] && echo "## ERROR: Please provide taxon ID(s) with -t flag" && exit 1
[[ "$kraken_main"  = "" ]] && echo "## ERROR: Please provide Kraken output file with -I flag" && exit 1
[[ ! -f $R1_in ]] && echo "## ERROR: Input file R1_in ($R1_in) does not exist" && exit 1
[[ ! -f $kraken_main ]] && echo "## ERROR: Kraken output file ($kraken_main) does not exist" && exit 1
[[ "$R1_in"  = "$R1_out" ]] && echo "## ERROR: R1 input and output filenames are the same: $R1_in" && exit 1

## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/krakentools-1.2

## Bash strict settings
set -euo pipefail

## Infer name of R2 input file
R2_in=${R1_in/_R1/_R2}

## Define output files
R1_out="$outdir"/$(basename "$R1_in")
R2_out="$outdir"/$(basename "$R2_in")

## Check input
[[ ! -f $R2_in ]] && echo "## ERROR: Input file R2_in ($R2_in) does not exist" && exit 1
[[ "$R2_in"  = "$R2_out" ]] && echo "## ERROR: R2 input and output filenames are the same: $R2_in" && exit 1

## Report
echo "## Taxonomic IDs to extract:          $taxids"
echo "## Kraken 'main' output file:         $kraken_main"
echo "## Input R1:                          $R1_in"
echo "## Input R2:                          $R2_in"
echo "## Output dir:                        $outdir"
echo "## Output R1:                         $R1_out"
echo "## Output R2:                         $R2_out"
echo -e "-----------------------------\n"

## Create output dir
mkdir -p "$outdir"


# EXTRACT READS ----------------------------------------------------------------
extract_kraken_reads.py \
    -t $taxids \
    -k "$kraken_main" \
    -s "$R1_in" -s2 "$R2_in" \
    -o "$R1_out" -o2 "$R2_out"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$R1_out" "$R2_out"
echo -e "\n## Done with script extract_kraken_reads.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo