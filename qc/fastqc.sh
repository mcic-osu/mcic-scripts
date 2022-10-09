#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=fastqc
#SBATCH --output=slurm-fastqc-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run FastQC for one or more FASTQ files."
    echo
    echo "Syntax: $0 -i <input-FASTQ> -o <output-dir> [fastq-file-1] [fastq-file-2] [...]"
    echo
    echo "Required options:"
    echo "    -i FILE    A single input FASTQ file (can be gzipped)"
    echi "               Alternatively, pass 1 or more FASTQ files as positional arguments at the end."
    echo "    -o DIR     Output directory"
    echo
    echo "Other options:"
    echo "    -h         Print this help message and exit"
    echo
    echo "Example command:"
    echo "    $0 -i data/fastq/my.fastq.gz -o results/fastqc"
    echo "    $0 -o results/fastqc data/fastq/sampleA.fastq.gz data/fastq/sampleB.fastq.gz"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
infiles=()
outdir=""

## Parse command-line options
while getopts ':i:o:h' flag; do
    case "${flag}" in
        i) infiles=("$OPTARG") ;;
        o) outdir="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done
shift "$(( OPTIND - 1 ))"

count=0
while [ "$*" != "" ]; do
    infiles[$count]=$1
    shift
    count=$(( count + 1 ))
done


# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/fastqc-0.11.9

## Input checks
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify output dir with -o"  >&2 && exit 1
[[ ! ${#infiles[@]} -gt 0 ]] && echo "## ERROR: Please specify one or more input files" >&2 && exit 1
for infile in ${infiles[*]}; do
    [[ ! -f "$infile" ]] && echo "## ERROR: Input file does not exist" >&2 && exit 1
done

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"/logs

## Report
echo "## Starting script fastqc.sh..."
date
echo
echo "## Input FASTQ file(s):    ${infiles[*]}"
echo "## Output dir:             $outdir"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
fastqc --outdir="$outdir" ${infiles[*]}


# WRAP UP ----------------------------------------------------------------------
sample_id=$(basename "$infile" | sed -E 's/.f?a?s?t?q.*//')
echo -e "\n## Listing output files:"
ls -lh "$outdir"/"$sample_id"*fastqc*
echo -e "\n## Done with script fastqc.sh"
date
echo
