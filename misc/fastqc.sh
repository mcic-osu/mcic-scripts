#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-fastqc-%j.out

set -u -e -o pipefail  # Bash strict settings

# Load any needed software
module load fastqc

# Processing command-line arguments
fastq_file="$1"
output_dir="$2"

# Create the output directory if it doesn't already exist
mkdir -p "$output_dir"

# Report:
date                              # Report date+time to time script
echo "Starting FastQC script..."  # Report what script is being run
echo "Running fastqc for file: $fastq_file"
echo "Output dir: $output_dir"
echo -e "---------\n\n"           # Separate from program output

# Run FastQC
fastqc --outdir="$output_dir" "$fastq_file"

# Report
echo -e "\n---------\nAll done!"  # Separate from program output
date                              # Report date+time to time script