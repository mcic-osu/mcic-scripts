#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --output=slurm-picrust-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8


# SETUP ------------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate picrust2-env-source

## Bash strict settings
set -euo pipefail

## Parse command-line args
biom=$1
fa=$2
outdir=$3

## Params for interactive testing
# biom=results/phyloseq/dadatax_filt.biom
# fa=results/ASV/main/ASVs.fa
# outdir=results/picrust
# n_cores=6

## Other parameters
n_proc=$SLURM_NTASKS

## Report
echo -e "\n## Starting script picrust.sh"
date
echo
echo "Input BIOM file:                    $biom"
echo "Input FASTA file:                   $fa"
echo "Output dir:                         $outdir"
echo "Number of processes:                $n_proc"
echo

## Test input
[[ ! -f $biom ]] && echo "## ERROR: Input BIOM file $biom does not exist" && exit 1
[[ ! -f $fa ]] && echo "## ERROR: Input BIOM file $biom does not exist" && exit 1


# RUN PICRUST ------------------------------------------------------------------
#? https://github.com/picrust/picrust2/wiki/Full-pipeline-script

picrust2_pipeline.py -s "$fa" \
                     -i "$biom" \
                     -o "$outdir" \
                     -p "$n_proc" \
                     --in_traits EC,KO \
                     --remove_intermediate \
                     --verbose

#? `--in_traits` - Comma-delimited list (with no spaces) of which gene families
#? to predict from this set: COG, EC, KO, PFAM, TIGRFAM. (default: EC,KO).


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir"

echo -e "\n## Done with script picrust.sh"
date
