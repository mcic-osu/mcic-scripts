#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-krona-%j.out

## Software
source ~/.bashrc
source activate krona-env

## Bash strict settings
set -euo pipefail

## Command-line args
infile=$1
outfile=$2

## Report
echo "## Starting script krona.sh..."
date
echo "## Command-line args:"
echo "## Input file: $infile"
echo "## Output file: $outfile"
echo -e "---------------------------------\n\n"

## Make outdir only if outfile is not in current dir (contains a "/")
if echo "$outfile" | grep -q "/"; then
    outdir=$(dirname "$outfile")
    echo "## Creating output dir $outdir"
    mkdir -p "$outdir"
fi

## Run Krona
ktImportTaxonomy -q 2 -t 3 "$infile" -o "$outfile"

## Report
echo "## Done with script krona.sh"
date
