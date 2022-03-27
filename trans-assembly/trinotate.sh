#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=trinotate
#SBATCH --output=slurm-trinotate-%j.out

## Report
echo -e "\n## Starting script trinotate.sh"
date
echo

## Load software
conda activate /users/PAS0471/jelmer/miniconda3/envs/trinotate-env
conda activate --stack /users/PAS0471/jelmer/miniconda3/envs/trinity-env # For script get_Trinity_gene_to_trans_map.pl
#> Note: it may be necessary to forego also loading the Trinity env in the later stages of the pipeline 

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in_org="$1"
config="$2"
outdir="$3"

## Process parameters - output files
GENE_TRANS_MAP=$outdir/GENE_TRANS_MAP

## Other paramaters
N_CORES=$SLURM_CPUS_PER_TASK

## Report
echo "## All args:                    $*"
echo
echo "## Input FASTA file:        $fa_in_org"
echo "## Input config file:       $config"
echo "## Output dir:              $outdir"
echo
echo "## Gene-to-transcript map:  $GENE_TRANS_MAP"
echo "## Number of cores:         $N_CORES"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"

## Copy assembly FASTA file to output dir
fa_in="$outdir"/$(basename "$fa_in_org")
if [ ! -f "$fa_in" ]; then
    echo "## Copying input FASTA to outdir..."
    cp -v "$fa_in_org" "$fa_in"
fi


# BUILD TRINOTATE DATABASE -----------------------------------------------------
if [ ! -f "$outdir"/Trinotate.sqlite ]; then
    echo "## Building Trinotate database..."
    cd "$outdir" || exit 1
    
    ## Build Trinotate DB:
    Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
    
    ## Uncompress and prepare the Pfam database for use with 'hmmscan':
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm

    ## Prepare the protein database for blast searches:
    makeblastdb -in uniprot_sprot.pep -dbtype prot
    
    cd - || exit 1
fi


# CREATE GENE-TO-TRANSCRIPTOME MAP ---------------------------------------------
if [ ! -f "$GENE_TRANS_MAP" ]; then
    echo -e "\n## Creating gene-to-transcript map..."
    get_Trinity_gene_to_trans_map.pl "$fa_in" > "$GENE_TRANS_MAP"
fi


# RUN TRINOTATE ----------------------------------------------------------------
echo -e "\n## Running Trinotate"

cd "$outdir" || exit 1

autoTrinotate.pl \
    --Trinotate_sqlite "$outdir"/Trinotate.sqlite \
    --transcripts "$fa_in" \
    --gene_to_trans_map "$GENE_TRANS_MAP" \
    --conf "$config" \
    --CPU "$N_CORES"

cd - || exit 1

# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script trinotate.sh"
date
echo

#? Info
# https://github.com/Trinotate/Trinotate
# https://github.com/Trinotate/Trinotate.github.io/wiki/Automated-Execution-of-Trinotate:-Running-computes-and-loading-results
# https://rnabio.org/module-07-trinotate/0007/02/01/Trinotate/