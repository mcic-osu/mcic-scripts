#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --job-name=fqconcat
#SBATCH --out=slurm-fqconcat-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Concatenate FASTQ files from the same sample but different lanes."
  echo
  echo "Syntax: $0 -i <input-dir> -o <output-dir>"
  echo
  echo "Required options:"
  echo "    -i DIR         Input dir containing FASTQ files, possibly in separate underlying dirs"
  echo "    -o DIR         Output dir for concatenated FASTQ files; all will be placed directly in this dir"
  echo
  echo "Other options:"
  echo "    -n INTEGER     Number of files expected per sample per read direction       [default: 2]"
  echo "    -p STRING      Globbing pattern for underlying dirs                         [default: '*' i.e. any dir]"
  echo "    -r STRING      Literal string or regex pattern to remove from file names    [default: none]"
  echo "    -h             Print this help message"
  echo
  echo "Example command:"
  echo "    $0 -i data/fastq/original -o data/fastq/concat -n 2"
  echo
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
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
  \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Input checks
[[ $indir = "" ]] && echo "## ERROR: Please specify input dir with -i" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify output dir with -o" && exit 1
[[ "$indir" = "$outdir" ]] && echo "## ERROR: Input dir $indir should not be the same as the output dir" && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir $indir does not exist" && exit 1

## Create output dir, if necessary
mkdir -p "$outdir"

## Report
echo "## Starting script fqconcat.sh"
date
echo
echo "## Input dir:                                              $indir"
echo "## Output dir:                                             $outdir"
[[ "$dir_pattern" != "" ]] && echo "## Dir globbing pattern:                                   $dir_pattern"
[[ "$remove_pattern" != "" ]] && echo "## Pattern/string to remove from file names:               $remove_pattern"
echo "## Nr files expected for each sample and read direction:   $nfiles"
echo -e "-------------------\n"


# CONCAT FASTQ FILES -----------------------------------------------------------
## Find per-sample directories
dirs=( $(find $indir -mindepth 1 -type d -name "*$dir_pattern*") )
echo "## Nr directories (samples) found:                         ${#dirs[@]}"
[[ ${#dirs[@]} -eq 0 ]] && echo "ERROR: No matching dirs found!" && exit 1
echo

## Loop over per-sample directories
for fdir in "${dirs[@]}"; do
    
    ## Count input files - exit if there are not 2 files per read direction
    n_R1=$(find "$fdir" -maxdepth 1 -name "*_R1_001.fastq.gz" | wc -l)
    n_R2=$(find "$fdir" -maxdepth 1 -name "*_R2_001.fastq.gz" | wc -l)

    ## Ignore dirs with no FASTQ files -- e.g. higher level dirs (!)
    if [ "$n_R1" != 0 ]; then
        echo -e "\n## Directory: $fdir"

        ## Check input
        [[ "$n_R1" != "$nfiles" ]] && echo "## ERROR: Did not find $nfiles R1 input files" && exit 1
        [[ "$n_R2" != "$nfiles" ]] && echo "## ERROR: Did not find $nfiles R2 input files" && exit 1

        ## Each dir should contain FASTQ files for a single sample.
        ## Extract the sample ID (e.g. "S9") from the 1st FASTQ file name.
        sample_id=$(basename "$(find "$fdir" -maxdepth 1 -name "*fastq.gz" | head -n 1)" .fastq.gz | sed -E 's/_L0[0-9][0-9]_.*//')
        sample_id="${sample_id/$remove_pattern/}"
        echo "## Sample ID: $sample_id"

        ## Make input dir read-only
        chmod a-w "$fdir"/*fastq.gz

        ## List input files
        echo "## R1 input files:"
        find "$fdir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh
        echo "## R2 input files:"
        find "$fdir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 du -sh
        
        ## Output files
        R1_out="$outdir"/"$sample_id"_R1_001.fastq.gz
        R2_out="$outdir"/"$sample_id"_R2_001.fastq.gz

        ## Concatenate FASTQ files
        find "$fdir" -name "*_R1_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$R1_out"
        find "$fdir" -name "*_R2_001.fastq.gz" -print0 | sort -z | xargs -0 -I{} cat {} > "$R2_out"

        ## List output files
        echo "## Output files:"
        du -h "$R1_out" "$R2_out"
    fi
    
    echo -e "------------------\n"
done

## Make concatenated files read-only to protect against accidental removal and overwriting
echo "## Making concatenated files read-only..."
chmod a-w "$outdir"/*fastq.gz

## Report
echo
echo "## Done with script fqconcat.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
