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
    echo "Syntax: $0 -i <FASTQ-in> -o <FASTQ-outdir> ..."
    echo
    echo "Required options:"
    echo "    -i FILE        (R1) FASTQ input file (if paired-end, R2 file name will be inferred)"
    echo "    -o DIR         Output dir"
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
while getopts ':i:o:q:l:sh' flag; do
    case "${flag}" in
        i) R1_in="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
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
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an outdir with -o" >&2 && exit 1
[[ "$R1_in" = "" ]] && echo "## ERROR: Please specify an input R1 FASTQ file with -i" >&2 && exit 1
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input R1 FASTQ file $R1_in does not exist" >&2 && exit 1

## Get number of threads
n_threads="$SLURM_CPUS_PER_TASK"

## Output dir for TrimGalore logs
outdir_trim="$outdir"/trimmed
outdir_fastqc="$outdir"/fastqc
outdir_logs="$outdir"/logs

## Get sample/file ID
extension=$(echo "$R1_in" | sed -E 's/.*(\.fa?s?t?q\.gz$)/\1/')

## Get R2 file, create input argument, define output files
if [ "$single_end" != "true" ]; then
    # Paired-end sequences
    R1_suffix=$(echo "$R1_in" | sed -E "s/.*(_R?1)_?[[:digit:]]*$extension/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    input_arg="--paired $R1_in $R2_in"
    [[ ! -f "$R2_in" ]] && echo "## ERROR: Input R2 FASTQ file $R2_in does not exist" >&2 && exit 1
    [[ "$R1_in" = "$R2_in" ]] && echo "## ERROR: Input R1 and R2 FASTQ files are the same file" >&2 && exit 1

    sample_id=$(basename "$R1_in" | sed -E "s/${R1_suffix}_?[[:digit:]]*${extension}//")
    R1_out="$outdir_trim"/"$sample_id"_R1.fastq.gz
    R2_out="$outdir_trim"/"$sample_id"_R2.fastq.gz
else
    # Single-end sequences
    input_arg="$R1_in"
    R1_suffix=""

    sample_id=$(basename "$R1_in" | sed "s/${R1_suffix}${extension}//")
    R1_out="$outdir_trim"/"$sample_id".fastq.gz
fi

## Report
echo -e "\n## Starting script trimgalore.sh"
date
echo
echo "## R1 input file:                     $R1_in"
echo "## Base output dir:                   $outdir"
echo
[[ "$single_end" != "true" ]] && echo "## R2 input file:                     $R2_in"
echo "## Sequence quality threshold:        $qual"
echo "## Minimum sequence length:           $len"
echo "## Sequences are single-end:          $single_end"
echo "## Sample ID:                         $sample_id"
echo
echo "## Output dir - FastQC:               $outdir_fastqc"
echo "## R1 output file:                    $R1_out"
[[ "$single_end" != "true" ]] && echo "## R2 output file:                    $R2_out"
echo -e "---------------------------\n"


# MAIN -------------------------------------------------------------------------
## Make output dirs
mkdir -p "$outdir_trim" "$outdir_fastqc" "$outdir_logs"

## Run Trim-Galore
trim_galore \
    --quality "$qual" \
    --length "$len" \
    --gzip \
    -j "$n_threads" \
    --output_dir "$outdir_trim" \
    --fastqc --fastqc_args "-t $n_threads --outdir $outdir_fastqc" \
    $input_arg

## Move output files
echo -e "\n## Moving output files..."
mv -v "$outdir_trim"/"$sample_id"*_trimming_report.txt "$outdir_logs"

if [ "$single_end" != "true" ]; then
    mv -v "$outdir_trim"/"$sample_id"*_val_1.fq.gz "$R1_out"
    mv -v "$outdir_trim"/"$sample_id"*_val_2.fq.gz "$R2_out"
else
    mv -v "$outdir_trim"/"$sample_id"*_trimmed.fq.gz "$R1_out"
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing FASTQ output files:"
ls -lh "$R1_out"
[[ "$single_end" != "true" ]] && ls -lh "$R2_out"
echo
echo -e "## Done with script trimgalore.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
