#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=mummer
#SBATCH --output=slurm-mummer-%j.out

# Mummer help
# https://mummer4.github.io/manual/manual.html

# Visualize Mummer output
# https://github.com/marianattestad/dot
# https://github.com/MariaNattestad/ribbon

# Load the Conda environment
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/mummer4

# Parse arguments
reference=$1
assembly=$2
outdir=$3

# Make paths absolute
[[ ! "$assembly" =~ ^/ ]] && assembly=$(realpath "$assembly")
[[ ! "$reference" =~ ^/ ]] && reference=$(realpath "$reference")

# Report
echo "Starting script dnadiff.sh"
date
echo "Input assembly:      $assembly"
echo "Input reference:     $reference"
echo "Output dir:          $outdir"
echo

# Move into the outdir
cd "$outdir" || exit 1

# Run dnadiff
dnadiff "$reference" "$assembly"

# Report
echo "Done with script"
date
