#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=effectorp
#SBATCH --output=slurm-effectorp-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run EffectorP to predict effectors in fungi/oomycetes."
  echo
  echo "Syntax: $0 -i <input-FASTA> -o <output-dir>..."
  echo
  echo "Required options:"
  echo      "-i STRING         Input protein FASTA file"
  echo      "-o STRING         Output directory"
  echo
  echo "Other options:"
  echo      "-a STRING         Other argument(s) to pass to EffectorP"
  echo      "-h                Print this help message"
  echo
  echo "Example:               $0 -i ref/aa.fa -o results/effectorp"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "EffectorP documentation: https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md"
  echo "EffectorP paper: https://apsjournals.apsnet.org/doi/10.1094/MPMI-08-21-0201-R"
  echo
}

## Option defaults
infile=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Check input
[[ "$infile" = "" ]] && echo "## ERROR: Please specify an input FASTA file with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output directory with -o" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input FASTA file (-i) $infile does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
EFFECTOR_P=/fs/project/PAS0471/jelmer/assist/2021-12_linda/software/EffectorP-3.0/EffectorP.py

## Bash strict mode
set -euo pipefail

## Make output dir
mkdir -p "$outdir"

## Define output files
table_out="$outdir"/table.txt
effector_fa="$outdir"/effectors.fa
noneffector_fa="$outdir"/noneffectors.fa

## Report
echo
echo "## Starting script effectorp.sh"
date
echo
echo "## Input AA FASTA:                        $infile"
echo "## Output dir:                            $outdir"
[[ "$more_args" != "" ]] && echo "## Other arguments to pass to EffectorP:    $more_args"
echo -e "--------------------\n"


# RUN LIFTOFF ------------------------------------------------------------------
echo "## Now running EffectorP..."
python "$EFFECTOR_P" \
    -i "$infile" \
    -o "$table_out" \
    -E "$effector_fa" \
    -N "$noneffector_fa" $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script effectorp.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
