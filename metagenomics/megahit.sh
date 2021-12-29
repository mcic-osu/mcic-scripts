#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-megahit-%j.out


# SETUP ------------------------------------------------------------------------
## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate megahit-env

## Bash strict mode
set -euo pipefail

## Command-line args
R1_in="$1"
outdir="$2"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                                       # Retrieve number of cores
mem=$(expr $(expr "$SLURM_MEM_PER_NODE" \* 1000000) - 1000000000)   # Retrieve memory

R2_in=${R1_in/_R1/_R2}
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?[0-9]).*fastq.gz/\1/')
R1_basename=$(basename "$R1_in" .fastq.gz)
sample_id=${R1_basename/"$R1_suffix"/}

contigs_fa="$outdir"/"$sample_id".contigs.fa      # This will be the output file that Megahit produces, containing all contigs 
single_contig_dir="$outdir"/single_contig_fa      # Dir with single-contig FASTA files

## Remove output dir if it exists (or Megahit will complain), make dir one level higher if needed
[[ -d "$outdir" ]] && echo "Removing $outdir" && rmdir "$outdir"
outdir_oneup=$(dirname "$outdir")
mkdir -p "$outdir_oneup"

## Report
echo "## Starting script megahit.sh..."
date
echo "## Command-line args:"
echo "## Input FASTQ file - R1:             $R1_in"
echo "## Output dir:                        $outdir"
echo "## Other parameters:"
echo "## Input FASTQ file - R2:             $R2_in"
echo "## Sample ID (used as output prefix): $sample_id"
echo "## Number of cores:                   $n_cores"
echo "## Memory in bytes:                   $mem"
echo -e "---------------------------\n\n"


# RUN MEGAHIT ------------------------------------------------------------------
echo "## Starting Megahit run..."
megahit \
    -1 "$R1_in" -2 "$R2_in" \
    -o "$outdir" \
    --out-prefix "$sample_id" \
    -t "$n_cores" -m "$mem" \
    --k-min 31 --k-max 111 --k-step 20 \
    --min-contig-len 500

## Create single-contig FASTA files
awk -v outdir="$single_contig_dir" -F " " \
    '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=outdir"/"ID".fa"} {print >> F}' < "$contigs_fa"

# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing file in output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script megahit.sh"
date