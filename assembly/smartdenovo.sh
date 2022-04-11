#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --output=slurm-smartdenovo-%j.out

## Load software
module load python
source activate /fs/project/PAS0471/jelmer/conda/smartdenovo-env

## Bash strict settings
set -euo pipefail

## Process command-line args
input=$1
prefix=$2

## Report
echo 
echo  "## Starting script smartdenovo.sh"
date
echo
echo "## Input:             $input"
echo "## Output prefix:     $prefix"
echo

## Create output dir
outdir=$(dirname "$prefix")
mkdir -p "$outdir"

## Generate a Makefile for smartdenovo to run
smartdenovo.pl \
    -p "$prefix" \
    -t "$SLURM_CPUS_PER_TASK" \
    -c 1 \
    "$input" \
    > "$prefix".mak

## Run SmartDenovo by running the Makefile
make -f "$prefix".mak

## Report
echo
echo "## Done with script smartdenovo.sh"
date

## Docs
# https://github.com/ruanjue/smartdenovo
# "After assembly, the raw unitigs are reported in file prefix.lay.utg and consensus unitigs in prefix.cns"
#> -c 1 => make consensus
#! -J min read length -- default = 5,000