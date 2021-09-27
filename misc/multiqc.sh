#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --out=slurm-multiqc-%j.out

# Software 
source ~/.bashrc
conda activate multiqc-env

# Bash strict mode
set -euo pipefail

# Command-line args
input_dir="$1"
output_dir=${2:-none} 

# Other variables
[[ "$output_dir" = "none" ]] && output_dir=$input_dir

# If necessary, create the output dir:
mkdir -p "$output_dir"

# Report before starting the pogram:
echo "Starting MultiQC script..."
date
echo "Running MultiQC for dir: $input_dir"
echo "Output dir: $output_dir"
echo -e "------------------\n\n"

# Run MultiQC
## --interactive will ensure interactive plots, regardless of number of samples
## --force will overwrite any old report
multiqc --interactive --force "$input_dir" -o "$output_dir"

# Report:
echo -e "\n Done with script multiqc.sh"
date

