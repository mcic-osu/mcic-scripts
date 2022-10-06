#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm-trimgalore-%j.out


# SETUP ------------------------------------------------------------------------
## Help
Help() {
    echo
    echo "$0: Run TrimGalore for one or a pair of FASTQ files"
    echo
    echo "Syntax: $0 -i <FASTQ-in> -o <FASTQ-outdir> -O <FastQC-outdir> ..."
    echo
    echo "Required options:"
    echo "    -i FILE        (R1) FASTQ input file (if paired-end, R2 file name will be inferred)"
    echo "    -o DIR         Trimmed FASTQ output dir"
    echo "    -O DIR         FastQC results output dir"
    echo
    echo "Other options:"
    echo "    -q INTEGER     Quality trimming threshold         [default: 20]"
    echo "    -l INTEGER     Minimum read length                [default: 20]"
    echo "    -s             Input is single-end                [default: paired-end]"
    echo "    -h             Print this help message and exit"
    echo
    echo "Example: $0 -i data/fastq/S01_R1.fastq.gz -o results/trimgalore -O results/fastqc/trimmed"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
qual=0                  # => 0 = no actual quality trimming
len=20                  # => 20 is also the TrimGalore default
single_end=false        # => paired-end by default

# Get command-line options:
while getopts ':i:o:O:q:l:sh' flag; do
    case "${flag}" in
        i) R1_in="$OPTARG" ;;
        o) outdir_trim="$OPTARG" ;;
        O) outdir_fastqc="$OPTARG" ;;
        q) qual="$OPTARG" ;;
        l) len="$OPTARG" ;;
        s) single_end=true ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/trimgalore-env

## Bash strict settings
set -euo pipefail

## Check input
[[ "$outdir_trim" = "" ]] && echo "## ERROR: Please specify an outdir for trimmed FASTQs with -o" >&2 && exit 1
[[ "$outdir_fastqc" = "" ]] && echo "## ERROR: Please specify an outdir for FastQC with -O" >&2 && exit 1
[[ "$R1_in" = "" ]] && echo "## ERROR: Please specify an input R1 FASTQ file with -i" >&2 && exit 1
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input R1 FASTQ file $R1_in does not exist" >&2 && exit 1

## Get number of threads
n_threads="$SLURM_CPUS_PER_TASK"

## Output dir for TrimGalore logs
logdir="$outdir_trim"/logs

## Get R2 file and create input argument
if [ "$single_end" != "true" ]; then
    R2_in=${R1_in/_R1_/_R2_}
    R2_id=$(basename "$R2_in" .fastq.gz)
    input_arg="--paired $R1_in $R2_in"
    [[ ! -f "$R2_in" ]] && echo "## ERROR: Input R2 FASTQ file $R2_in does not exist" >&2 && exit 1
else
    input_arg="$R1_in"
fi

## Get sample/file ID
R1_id=$(basename "$R1_in" .fastq.gz)

## Make output dirs
mkdir -p "$outdir_trim" "$logdir" "$outdir_fastqc"

## Report
echo -e "\n## Starting script trimgalore.sh"
date
echo
echo "## R1 input file:                 $R1_in"
echo "## Output dir - trimmed FASTQs:   $outdir_trim"
echo "## Output dir - FastQC:           $outdir_fastqc"
echo "## Sequence quality threshold:    $qual"
echo "## Minimum sequence length:       $len"
[[ "$single_end" != "true" ]] && echo "## R2 input file:                 $R2_in"
echo -e "---------------------------\n\n"


# RUN TRIMGALORE AND PROCESS OUTPUT --------------------------------------------
## Run Trim-Galore
trim_galore \
    --quality "$qual" \
    --length "$len" \
    --gzip \
    -j "$n_threads" \
    --output_dir "$outdir_trim" \
    --fastqc --fastqc_args "-t $n_threads --outdir $outdir_fastqc" \
    $input_arg

## Move log files
mv "$outdir_trim"/*_trimmming_report.txt "$logdir"

## Move FASTQ files
if [ "$single_end" != "true" ]; then
    mv "$outdir_trim"/"$R1_id"_val_1.fq.gz "$outdir_trim"/"$R1_id".fastq.gz
    mv "$outdir_trim"/"$R2_id"_val_2.fq.gz "$outdir_trim"/"$R2_id".fastq.gz
else
    mv "$outdir_trim"/"$R1_id"_trimmed.fq.gz "$outdir_trim"/"$R1_id".fastq.gz
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir_trim"/"$R1_id".fastq.gz "$outdir_trim"/"$R2_id".fastq.gz
echo -e "\n## Done with script trimgalore.sh"
date
echo
