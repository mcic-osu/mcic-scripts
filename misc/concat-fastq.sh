#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --out=slurm-concat-fastq-%j.out

## Bash strict settings
set -euo pipefail

## Help
Help() {
  echo
  echo "## $0: Concatenate FASTQ files from different lanes."
  echo
  echo "## Syntax: $0 -i <input-dir> -o <output-dir> [-p <dir-glob-pattern>] [-r <remove-pattern>] [-n <n-files> ] [-h]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   Input directory containing FASTQ files, possibly in separate underlying dirs (REQUIRED)"
  echo "## -o STR   Output directory for concatenated FASTQ files; all will be placed directly in this dir (REQUIRED)"
  echo "## -p STR   Globbing pattern for underlying dirs (default: '*' i.e. any dir)"
  echo "## -r STR   Literal string or regex pattern to remove from file names (default: none)"
  echo "## -n INT   Number of files expected per sample per read direction (default: 2)"
  echo "## Example: $0 -i data/fastq/original -o data/fastq/concat -n 2"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
indir=""
outdir=""
dir_pattern="*"
remove_pattern=""
nfiles=2

## Parse command-line options
while getopts ':i:o:p:r:n:h' flag; do
  case "${flag}" in
  i) indir="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  p) dir_pattern="$OPTARG" ;;
  r) remove_pattern="$OPTARG" ;;
  n) nfiles="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Input checks
[[ $indir = "" ]] && echo "## ERROR: Please specify input dir with -i" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify input dir with -i" && exit 1
[[ "$indir" = "$outdir" ]] && echo "## ERROR: Input dir should not be the same as the output dir" && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir does not exist" && exit 1

## Create output dir, if necessary
mkdir -p "$outdir"

## Report
echo "## Starting script concat-fastq.sh"
date
echo "## Input dir:: $indir"
echo "## Output dir: $outdir"
echo "## Dir globbing pattern: $dir_pattern"
echo "## Pattern/string to remove from file names: $remove_pattern"
echo "## Number of files expected for each sample and read direction: $nfiles"
echo -e "-------------------\n\n"

## Find per-sample directories
dirs=( $(find $indir -mindepth 1 -type d -name "$dir_pattern") )
echo "Number of directories (samples): ${#dirs[@]}"

## Loop over per-sample directories
for fdir in "${dirs[@]}"; do
    echo -e "\n## Directory: $fdir"
    
    ## Each dir should contain FASTQ files for a single sample.
    ## First, extract the sample ID (e.g. "S9") from the 1st FASTQ file name.
    ## Then, zero-pad single-digit numbers (e.g. "S9" -> "S09") to ensure proper ordering.
    sample_id=$(basename "$(find "$fdir" -name "*fastq.gz" | head -n 1)" .fastq.gz | sed -E -e 's/_L0[0-9][0-9]_.*//')
    sample_id="${sample_id/$remove_pattern/}"
    echo "## Sample ID: $sample_id"

    ## Count input files - exit if there are not 2 files per read direction
    n_R1=$(find "$fdir" -name "*_R1_001.fastq.gz" | wc -l)
    n_R2=$(find "$fdir" -name "*_R2_001.fastq.gz" | wc -l)
    [[ "$n_R1" != "$nfiles" ]] && echo "## ERROR: Did not find two R1 input files" && exit 1
    [[ "$n_R2" != "$nfiles" ]] && echo "## ERROR: Did not find two R2 input files" && exit 1

    ## List input files
    echo "## R1 input files:"
    find "$fdir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh
    echo "## R2 input files:"
    find "$fdir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh
    
    ## Output files
    R1_out="$outdir"/"$sample_id"_R1_001.fastq.gz
    R2_out="$outdir"/"$sample_id"_R2_001.fastq.gz
    echo "## Output files:"
    echo "$R1_out"
    echo "$R2_out"

    ## Concatenate FASTQ files
    find "$fdir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$R1_out"
    find "$fdir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$R2_out"

    ## List output files
    echo "## Output files:"
    du -h "$outdir"/"$sample_id"_R1.fastq.gz
    du -h "$outdir"/"$sample_id"_R2.fastq.gz

    echo -e "---------\n"
done

## Make concatenated files read-only to protect against accidental removal and overwriting
echo "## Making concatenated files read-only..."
chmod a-w "$outdir"/*fastq.gz

## Report
echo -e "\n## Done with script concat-fastq.sh."
date