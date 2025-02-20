#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --job-name=fq2fa
#SBATCH --output=slurm-fq2fa-%j.out

# Convert FASTQ to FASTA with seqtk

# Command-line args
fq_in=$1
fa_out=$2

# Software
module load miniconda3
source activate /fs/ess/PAS0471/jelmer/conda/seqtk

# Bash strict settings
set -euo pipefail

# Report
echo "# Starting script fq2fa.sh"
echo "# Input FASTQ file:       $fq_in"
echo "# Output FASTA file:      $fa_out"
date

# Convert FASTQ to FASTA
seqtk seq -a "$fq_in" > "$fa_out"

# Report
echo "# Done with script fq2fa.sh"
date
echo
