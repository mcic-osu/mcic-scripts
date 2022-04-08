#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=sortmerna
#SBATCH --output=slurm-sortmerna-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Run SortMeRNA to sort FASTQ into rRNA-derived and other reads"
    echo
    echo "## Syntax: $0 -i <R1-FASTQ-file> -o <output-dir> ..."
    echo 
    echo "## Required options:"
    echo "## -i STRING       Input R1 FASTQ file"
    echo "## -o STRING       Output directory"
    echo
    echo "## Other options:"
    echo "## -h              Print this help message and exit"
    echo
    echo "## Example: $0 -i data/A1_R1_001.fastq.gz -o results/sortmerna"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
R1=""
outdir=""

## Get command-line options
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) R1="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Report
echo "## Starting script sortmerna.sh"
date
echo

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/sortmerna-env
BBTOOLS_ENV=/users/PAS0471/jelmer/miniconda3/envs/bbmap-env # Used at end

## Bash strict settings
set -euo pipefail

## Infer the name of the R2 file
R2=${R1/_R1_/_R2_}

## Infer the sampleID and define the fill output dir
sampleID=$(basename "$R1" | sed 's/_R1.*//')
outdir_full="$outdir"/"$sampleID"

repo_dir="$outdir"/sortmerna_repo

out_aligned="$outdir"/aligned_tmp/"$sampleID"
out_nonaligned="$outdir"/nonaligned_tmp/"$sampleID"

R1_aligned="$outdir"/aligned/"$sampleID"_R1_001.fastq.gz
R2_aligned="$outdir"/aligned/"$sampleID"_R2_001.fastq.gz
R1_nonaligned="$outdir"/nonaligned/"$sampleID"_R1_001.fastq.gz
R2_nonaligned="$outdir"/nonaligned/"$sampleID"_R2_001.fastq.gz

## Reference FASTA files (to be downloaded)
ref_18s="$repo_dir"/data/rRNA_databases/silva-euk-18s-id95.fasta
ref_28s="$repo_dir"/data/rRNA_databases/silva-euk-28s-id98.fasta

## Report
echo "## R1 FASTQ file:              $R1"
echo "## R2 FASTQ file:              $R2"
echo "## Output dir:                 $outdir_full"
echo "## SortMeRNA repo dir:         $repo_dir"
echo "## 18S reference file:         $ref_18s"
echo "## 28S reference file:         $ref_28s"         
echo -e "---------------------------\n"

## Check input
[[ ! -f "$R1" ]] && echo "## ERROR: R1 FASTQ file $R1 not found" >&2 && exit 1
[[ ! -f "$R2" ]] && echo "## ERROR: R2 FASTQ file $R2 not found" >&2 && exit 1

## Make output dir if needed
mkdir -p "$outdir"/aligned_tmp "$outdir"/nonaligned_tmp "$outdir"/aligned "$outdir"/nonaligned


# GET DATABASE FILES -----------------------------------------------------------
## Clone sortmerna repo to get db FASTA files
if [ ! -d "$repo_dir" ]; then
    mkdir -p "$repo_dir"
    echo "## Cloning sortmerna repo..."
    git clone https://github.com/biocore/sortmerna "$repo_dir"
fi

## Check that db files are there
[[ ! -f "$ref_18s" ]] && echo "## ERROR: 18s reference FASTA file $ref_18s not found" >&2 && exit 1
[[ ! -f "$ref_28s" ]] && echo "## ERROR: 28s reference FASTA file $ref_28s not found" >&2 && exit 1


# RUN SortMeRNA ----------------------------------------------------------------
echo -e "## Starting SortMeRNA run....\n"
sortmerna \
    --ref "$ref_18s" \
    --ref "$ref_28s" \
    --reads "$R1" \
    --reads "$R2" \
    --fastx \
    --aligned "$out_aligned" \
    --other "$out_nonaligned" \
    --workdir "$outdir_full" \
    --threads "$SLURM_CPUS_PER_TASK"


# CONVERTING INTERLEAVED FASTQ BACK TO SEPARATED -------------------------------
set +u
conda deactivate
source activate "$BBTOOLS_ENV"
set -u

echo -e "\n## Converting R1 back to paired files..."
reformat.sh \
    in="$out_aligned".fq.gz \
    out1="$R1_aligned" \
    out2="$R2_aligned"

echo -e "\n## Converting R2 back to paired files..."
reformat.sh \
    in="$out_nonaligned".fq.gz \
    out1="$R1_nonaligned" \
    out2="$R2_nonaligned"

[[ -f "$R1_aligned" ]] && rm "$out_aligned".fq.gz
[[ -f "$R1_nonaligned" ]] && rm "$out_nonaligned".fq.gz


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$R1_aligned" "$R2_aligned" "$R1_nonaligned" "$R2_nonaligned"
echo -e "\n## Done with script sortmerna.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo