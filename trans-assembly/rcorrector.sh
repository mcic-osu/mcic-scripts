#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=rcorrector
#SBATCH --output=slurm-rcorr-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Run rcorrector for a directory of FASTQ files."
  echo
  echo "## Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING     Input directory with FASTQ files"
  echo "## -o STRING     Output directory for corrected FASTQ FILES"
  echo
  echo "## Other options:"
  echo "## -s INTEGER    Start at specified step           [default: '0']"
  echo "                 Should be 0 unless restarting an interrupted/failed run"
  echo "## -h            Print this help message and exit"
  echo
  echo "## Example: $0 -i data/fastq/ -o results/rcorrector"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
indir=""
outdir=""
start_at_step=0

## Get command-line parameter values
while getopts ':i:o:s:h' flag; do
    case "${flag}" in
    i) indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    s) start_at_step="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script rcorrector.sh"
date
echo

## Test parameter values
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-i) $indir does not exist" >&2 && exit 1


# SOFTWARE ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/rcorrector-env


# OTHER SETUP ------------------------------------------------------------------
## Bash strict mode
set -euo pipefail

## Process parameters
R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')

## Report
echo "## Input dir:           $indir"
echo "## Output dir:          $outdir"
echo
echo "## List of R1 files:    $R1_list"
echo "## List of R2 files:    $R2_list"
echo
echo "## Start at step:       $start_at_step"
echo -e "--------------------------\n"

## Create output dirs if needed
mkdir -p "$outdir"


# RUN RCORRECTOR ---------------------------------------------------------------
run_rcorrector.pl \
    -t "$SLURM_CPUS_PER_TASK" \
    -od "$outdir" \
    -stage "$start_at_step" \
    -1 "$R1_list" -2 "$R2_list"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n--------------------------"
echo "## Listing file in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script rcorrector.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
