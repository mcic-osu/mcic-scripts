#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=snippy
#SBATCH --output=slurm-snippy-%j.out

# FUNCTIONS --------------------------------------------------------------------
## Help function
print_help() {
    echo
    echo "======================================================================"
    echo "$0: Run snippy to align FASTQ files to a genome and find SNPs"
    echo "======================================================================"
    echo "USAGE:"
    echo "------------------"
    echo "$0 -i <input-table> -r <reference-FASTA> -o <output-dir> ..."
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i FILE           Input TSV with paths to FASTQ files"
    echo "                      Should have 3 columns (no header): genome ID, forward reads FASTQ, reverse reads FASTQ"
    echo "    -r FILE           Reference genome FASTA file"
    echo "    -o DIR            Output dir"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -a STRING         Other argument(s) to pass to Snippy"
    echo "    -v                Print the version and exit"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "sbatch $0 -i snippy-input.tsv -r results/assembly/genome.fa -o results/snippy"
    echo
    echo "DOCUMENTATION:"
    echo "------------------"
    echo "Repo/documentation:   https://github.com/tseemann/snippy"
    echo
}

## Load software
load_software() {
    module load python/3.6-conda5.2
    source activate /fs/project/PAS0471/jelmer/conda/snippy-4.6.0
}

## Print the version
print_version() {
    load_software
    treetime --version
}


# PARSE OPTIONS ----------------------------------------------------------------
## Option defaults
input_table=""
ref_fasta=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:r:o:a:vh' flag; do
    case "${flag}" in
        i) input_table="$OPTARG" ;;
        r) ref_fasta="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        v) print_version && exit 0 ;;
        h) print_help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ "$input_table" = "" ]] && echo "## ERROR: Please specify an input table with -i" >&2 && exit 1
[[ "$ref_fasta" = "" ]] && echo "## ERROR: Please specify a reference FASTA file file with -r" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$input_table" ]] && echo "## ERROR: Input file $input_table does not exist" >&2 && exit 1
[[ ! -f "$ref_fasta" ]] && echo "## ERROR: Input file $ref_fasta does not exist" >&2 && exit 1

## Make path to reference absolute
[[ ! "$ref_fasta" =~ ^/ ]] && ref_fasta="$PWD"/"$ref_fasta"

## Report
echo "=========================================================================="
echo "## Starting script snippy.sh"
date
echo "=========================================================================="
echo
echo "## Input TSV with paths to FASTQ files:       $input_table"
echo "## Reference genome FASTA file:               $ref_fasta"
echo "## Output dir:                                $outdir"
[[ $more_args != "" ]] && echo "## Other arguments for Snippy:              $more_args"
echo
echo "## Listing the input files:"
ls -lh "$input_table" "$ref_fasta"
echo
echo "## Showing the contents of the input table file:"
cat -n "$input_table"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
## Create output dir
mkdir -p "$outdir"/logs

echo "## Generating Snippy commands with snippy-multi"
set -o xtrace
snippy-multi \
    "$input_table" \
    --ref "$ref_fasta" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    > "$outdir"/runme.sh
set +o xtrace
#? use `--prefix` option?

echo -e "\n## Showing contents of generated runme.sh:"
cat -n "$outdir"/runme.sh

echo -e "\n## Now running Snippy..."
cd "$outdir" || exit 1
bash runme.sh


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Snippy version used:"
snippy --version
echo
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo 
echo "## Done with script snippy.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
