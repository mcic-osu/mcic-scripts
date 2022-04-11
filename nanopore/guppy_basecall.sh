#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm-guppy_basecall-%j.out

# SETUP ------------------------------------------------------------------------
## Load software 
GUPPY_DIR=/fs/project/PAS0471/jelmer/software/guppy-6.0.1
GUPPY_BIN="$GUPPY_DIR"/bin/guppy_basecaller

## Bash strict settings
set -ueo pipefail

## Command-line args
indir=$1
outdir=$2
config=$3
barcode_kit=$4
min_qscore=$5

## Input checks
[[ "$#" != 5 ]] && echo "## ERROR: Please provide 5 arguments; you provided $#" && exit 1
#[[ ! -d "$indir" ]] && echo "## ERROR: Input dir does not exist" && exit 1

## Other parameters
N_CORES=$SLURM_CPUS_PER_TASK

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"

## Report
echo "## Starting script guppy_basecall.sh..."
date
echo
echo "## Input dir with FASTQ files:   $indir"
echo "## Output dir:                   $outdir"
echo "## ONT config name:              $config"
echo "## Barcode kit name:             $barcode_kit"
echo "## Minimum qual score to PASS:   $min_qscore"
echo "## Number of cores (threads):    $N_CORES"
echo -e "--------------------\n"


# RUN GUPPY -------------------------------------------------------------------
$GUPPY_BIN \
    --input_path "$indir" \
    --save_path "$outdir" \
    --config "$config" \
    --barcode_kits "$barcode_kit" \
    --min_qscore "$min_qscore" \
    --trim_barcodes \
    --compress_fastq \
    --records_per_fastq 0 \
    --cpu_threads_per_caller "$N_CORES" \
    --num_barcode_threads "$N_CORES" \
    --num_callers 1


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing files in output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script guppy_basecall.sh"
date


## To print config names for flowcell + kit combs:
#$ guppy_basecaller --print_workflows
## For kit LSK109 and flowcell FLO-MIN106, this is dna_r9.4.1_450bps_hac