#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=300
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-kraken-run-%j.out

## Help function
Help() {
    echo
    echo "## $0: Run Kraken2."
    echo
    echo "## Syntax: $0 -i <input-sequence-file> -o <output-dir> -d <kraken-db-dir> ..."
    echo 
    echo "## Required options:"
    echo "## -i     Input sequence file (FASTA, single-end FASTQ, or R1 from paired-end FASTQ)"
    echo "## -o     Output dir"
    echo "## -d     Kraken database dir"
    echo
    echo "## Other options:"
    echo "## -h     Print help."
    echo "## -n     Add taxonomic names to 'main' file (not compatible with Krona)"
    echo "## Example: $0 -i refdata/kraken/my-db -u https://genome.fa -g refdata/kraken/my-db/genomes"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env

## Bash strict settings
set -euo pipefail

## Option defaults
infile=""
outdir=""
krakendb_dir=""
add_names=false

## Get command-line options
while getopts 'i:o:d:nh' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) krakendb_dir="$OPTARG" ;;
    n) add_names=true ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process options
[[ "$infile" = "" ]] && echo "ERROR: must specify an input file with -i" && exit 1
[[ "$outdir" = "" ]] && echo "ERROR: must specify an output dir with -o" && exit 1
[[ "$krakendb_dir" = "" ]] && echo "ERROR: must specify an Kraken DB dir with -d" && exit 1

## Report
echo "## Starting script kraken-run..."
date
echo "## Command-line args:"
echo "## Input file: $infile"
echo "## Output dir: $outdir"
echo "## Kraken db dir: $krakendb_dir"
echo "## Add tax. names (Krona compatibility): $add_names"
echo

## Create output dir
mkdir -p "$outdir"

## Add tax. names or onot
if [ "$add_names" = true ]; then
    names_arg="--use-names "
    #? report-minimizer-data: see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information
else
    names_arg=""
fi

## Make sure input file argument is correct based on file type 
if [[ "$infile" =~ \.fastq.gz$ ]]; then
    R1_in="$infile"
    R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?[0-9]).*fastq.gz/\1/')
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    R1_basename=$(basename "$R1_in" .fastq.gz)
    sample_ID=${R1_basename/"$R1_suffix"/}

    if [[ -f $R2_in ]]; then
        echo "## Input is: paired FASTQ files"
        echo "## Input FASTQ file - R1: $R1_in"
        echo "## Input FASTQ file - R2: $R2_in"
        infile_arg="--gzip-compressed --paired $R1_in $R2_in"
    else
        echo -e "\n## Input is: single-end FASTQ file"
        infile_arg="--gzip-compressed $R1_in"
    fi
else
    echo -e "\n## Input is: FASTA file"
    infile_basename=$(basename "$infile")
    sample_ID=${infile_basename%%.*}
    infile_arg="$infile"
fi

## Define output files
outfile_main="$outdir"/"$sample_ID"_main.txt
outfile_report="$outdir"/"$sample_ID"_report.txt

## Report
echo "## add_names / add_names arg: $add_names $names_arg"
echo "## Input file arg: $infile_arg"
echo "## Sample ID: $sample_ID"
echo "## Output file - main: $outfile_main"
echo "## Output file - report: $outfile_report"
echo -e "------------------------------\n"


# RUN KRAKEN -------------------------------------------------------------------
echo "## Starting Kraken2 run..."

kraken2 ${names_arg}--threads "$SLURM_CPUS_ON_NODE" \
    --report-minimizer-data \
    --db "$krakendb_dir" \
    --report "$outfile_report" \
    ${infile_arg}>"$outfile_main"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outfile_main" "$outfile_report"
echo -e "\n## Done with script kraken-run.sh."
date
