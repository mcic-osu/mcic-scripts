#!/bin/bash

## Bash strict settings
set -euo pipefail

## Constants
SUBSAMPLE_SCRIPT="mcic-scripts/utils/fqsub.sh"

# PARSE ARGS -------------------------------------------------------------------
## Help function
Help()
{
   echo
   echo "fqsub_dir.sh: script to subsample a dir with FASTQ files using seqtk"
   echo
   echo "Syntax: fqsub_dir.sh -i <input-dir> -o <output-dir> ..."
   echo "NOTE: Don't submit this script to the SLURM queue -- it will itself spawn jobs instead."
   echo
   echo "Required options:"
   echo "   -i     Input dir"
   echo "   -o     Output dir"
   echo
   echo "Other options:"
   echo "   -s     Single-end sequences         [default: paired-end]"
   echo "   -n     Number of reads              [default: 100,000]"
   echo "   -p     Proportion of reads          [alternative to '-n', not applied by default]"
   echo "   -x     Sample ID pattern (select only matching filenames)"
   echo "   -h     Print this help message and exit"
   echo
}

## Option defaults
n_reads=100000
prop_reads=""
sample_pattern=""
single_end=false

## Parse command-line options
while getopts ':i:o:n:p:x:sh' flag; do
	case "${flag}" in
		i)  indir="$OPTARG" ;;
		o)  outdir="$OPTARG" ;;
		n)  n_reads="$OPTARG" ;;
		p)  prop_reads="$OPTARG" ;;
		s)  single_end=true ;;
		x)  sample_pattern="$OPTARG" ;;
		h)  Help && exit 0 ;;
		\?) echo "## ERROR: Invalid option" >&2 && exit 1 ;;
		:)  echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
	esac
done

## Check input
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir $indir does not exist" >&2 && exit 1
[[ "$indir" = "$outdir" ]] && echo "## ERROR: Input dir $indir cannot be the same as the output dir $outdir" >&2 && exit 1

## Report
echo
echo "## Starting script fqsub_dir.sh"
date
echo "## Input dir:                    $indir"
echo "## Output dir:                   $outdir"
[[ $n_reads != "" ]] && echo "## Number of reads to keep:      $n_reads"
[[ $prop_reads != "" ]] && echo "## Proportion of reads to keep:  $prop_reads"
echo "## Are reads single-end?         $single_end"
[[ $sample_pattern != "" ]] && echo "## Sample ID pattern:            $sample_pattern"
echo -e "------------------------\n"

## Make output dir if needed
mkdir -p "$outdir"


# RUN SCRIPT FOR EACH PAIR OF FASTQ FILES --------------------------------------
if [ "$single_end" = false ]; then

    for R1 in "$indir"/${sample_pattern}*_R1*q.gz; do
        R1=$(basename "$R1")
        R2=${R1/_R1/_R2}

        ## Report input files:
        echo "## R1 input file:"
        ls -lh "$indir/$R1"
        echo "## R2 input file:"
        ls -lh "$indir/$R2"

        ## SLURM log file
        sample_id=$(basename "$R1" .fastq.gz)
        log=slurm-fqsub-"$sample_id"-%j.out

        sbatch -o "$log" "$SUBSAMPLE_SCRIPT" \
            -i "$indir/$R1" -I "$indir/$R2" \
            -o "$outdir/$R1" -O "$outdir/$R2" \
            -n "$n_reads" -p "$prop_reads"
    
        echo -e "---------------------------\n"
    done

else

    for R1 in "$indir"/${sample_pattern}*q.gz; do
        R1=$(basename "$R1")

        ## Report input file:
        echo "## Input file:"
        ls -lh "$indir/$R1"

        ## SLURM log file
        sample_id=$(basename "$R1" .fastq.gz)
        log=slurm-fqsub-"$sample_id"-%j.out

        sbatch -o "$log" "$SUBSAMPLE_SCRIPT" \
            -i "$indir/$R1" \
            -o "$outdir/$R1" \
            -n "$n_reads" -p "$prop_reads"
    
        echo -e "---------------------------\n"
    done

fi

# REPORT AND FINALIZE --------------------------------------------------------
echo -e "## Done with script fqsub_dir.sh"
date
