#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-picrust-%j.out

# SETUP ------------------------------------------------------------------------
# Load software
module load miniconda3
conda activate XXX #TODO

# Bash strict settings
set -euo pipefail

# Parse command-line args
biom=$1             # .biom file
fa=$2               # FASTA file with ASVs
outdir=$3           # Output dir

# Report
echo -e "\n# Starting script picrust.sh"
date
echo
echo "Input BIOM file:                    $biom"
echo "Input FASTA file:                   $fa"
echo "Output dir:                         $outdir"
echo

# Test input
[[ ! -f $biom ]] && echo "# ERROR: Input BIOM file $biom does not exist" && exit 1
[[ ! -f $fa ]] && echo "# ERROR: Input FASTA file $fa does not exist" && exit 1

# RUN PICRUST ------------------------------------------------------------------
#? https://github.com/picrust/picrust2/wiki/Full-pipeline-script
#? `--in_traits` - Comma-delimited list (with no spaces) of which gene families
#? to predict from this set: COG, EC, KO, PFAM, TIGRFAM. (default: EC,KO).

picrust2_pipeline.py \
    -s "$fa" \
    -i "$biom" \
    -o "$outdir" \
    -p 8 \
    --in_traits EC,KO \
    --remove_intermediate \
    --verbose

# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n# Listing output files:"
ls -lh "$outdir"
echo -e "\n# Done with script picrust.sh"
date
