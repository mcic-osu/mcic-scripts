#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH --output=slurm-quast-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run QUAST to check the quality of a genome assembly."
  echo
  echo "Syntax: $0 [ -i <input-file> / -d <input-dir> ] -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Input file (genome assembly nucleotide FASTA)"
  echo "                          Specify the input with EITHER '-i' (single file) or '-d' (all files in a dir)"  
  echo "    -d STRING         Input dir"
  echo "                          Note: Only files with the extension '.fasta' will be included"
  echo "                          Specify the input with EITHER '-i' (single file) or '-d' (all files in a dir)"  
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -r STRING         Reference genome nucleotide FASTA file    [default: no reference genome]"
  echo "    -g STRING         Reference genome GFF file                 [default: no reference genome]"
  echo "    -1 STRING         FASTQ file with forward reads             [default: no reads]"
  echo "    -2 STRING         FASTQ file with reverse reads             [default: no reads]"
  echo "    -a STRING         Other argument(s) to pass to QUAST"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:          $0 -i results/assembly -o results/quast"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Quast documentation: https://github.com/ablab/quast"
  echo "Quast paper: https://academic.oup.com/bioinformatics/article/29/8/1072/228832"
  echo
}

## Option defaults
infile=""
indir=""
outdir=""
ref_fna=""
ref_gff=""
R1=""
R2=""
more_args=""

## Parse command-line options
while getopts '1:2:i:d:o:r:g:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    d) indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    r) ref_fna="$OPTARG" ;;
    g) ref_gff="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Check input
[[ "$indir" = "" && "$infile" = "" ]] && echo "## ERROR: Please specify either an input file with -i, or an input dir with -d" >&2 && exit 1
[[ "$indir" != "" && "$infile" != "" ]] && echo "## ERROR: Please specify either an input file with -i, or an input dir with -d, and not both!" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ "$R2" != "" && "$R1" = "" ]] && echo "## ERROR: When specifying a FASTQ file with reverse reads ('-2'), please also specify a FASTQ file with forward reads ('-1')" >&2 && exit 1
[[ "$infile" != "" && ! -f "$infile" ]] && echo "## ERROR: Input file (-i) $infile does not exist or is not a regular file" >&2 && exit 1
[[ "$indir" != "" && ! -d "$indir" ]] && echo "## ERROR: Input dir (-d) $indir does not exist or is not a directory" >&2 && exit 1

## Load the software
module load python
source activate /fs/project/PAS0471/jelmer/conda/quast-5.0.2

## Bash strict settings
set -euo pipefail

## Build argument for assembly input
if [[ $indir != "" ]]; then
    infile_arg="$indir/*.fasta"
else
    infile_arg="$infile"
fi

## If the input is a single file, make a separate output dir
if [[ $infile != "" ]]; then
    sampleID=$(basename "$infile" | sed -E 's/.fn?as?t?a?//')
    outdir=$outdir/"$sampleID"
fi

## Build argument for FASTQ files
if [[ $R1 != "" ]]; then
    R1_arg="-1 $R1"
else
    R1_arg=""
fi

if [[ $R2 != "" ]]; then
    R2_arg="-2 $R2"
else
    R2_arg=""
fi

## Build arguments for ref genome
if [[ $ref_fna != "" ]]; then
    ref_fna_arg="-r $ref_fna"
else
    ref_fna_arg=""
fi

if [[ $ref_gff != "" ]]; then
    ref_gff_arg="-g $ref_gff"
else
    ref_gff_arg=""
fi

## Make output dir
mkdir -p "$outdir"

## Report
echo
echo "## Starting script quast.sh"
date
echo
[[ $infile != "" ]] && echo "## Input file:                           $infile"
[[ $indir != "" ]] && echo "## Input dir:                            $indir"
echo "## Output dir:                           $outdir"
[[ $ref_fna != "" ]] && echo "## Reference FASTA file:                 $ref_fna"
[[ $ref_gff != "" ]] && echo "## Reference GFF file:                   $ref_gff"
[[ $R1 != "" ]] && echo "## R1 FASTQ file:                        $R1"
[[ $R2 != "" ]] && echo "## R2 FASTQ file:                        $R2"
[[ $more_args != "" ]] && echo "## Other arguments to pass to QUAST:     $more_args"
if [[ $indir != "" ]]; then
    echo -e "\n## Listing input FASTA files:"
    ls -lh "$indir"/*.fasta
fi
echo -e "--------------------\n"


# RUN QUAST --------------------------------------------------------------------
quast.py \
    --conserved-genes-finding \
    -o "$outdir" \
    -t "$SLURM_CPUS_PER_TASK" $ref_fna_arg $ref_gff_arg $R1_arg $R2_arg $more_args \
    $infile_arg


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script quast.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
