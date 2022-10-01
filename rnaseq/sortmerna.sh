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
    echo "## -r              Directory with SortMeRNA repo (for reference FASTA files)  [default: download repo]"
    echo "## -d              Don't 'de-interleave' output FASTQ file                    [default: de-interleave]"
    echo "## -h              Print this help message and exit"
    echo
    echo "## Example: $0 -i data/A1_R1_001.fastq.gz -o results/sortmerna"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "## Output:"
    echo "##   - Aligned sequences will be placed in <output-dir>/mapped"
    echo "##   - Non-aligned sequences will be placed in <output-dir>/unmapped"
    echo "## Output sequence files will keep sample identifiers,"
    echo "##   so the script can be run for multiple samples using the same <output-dir> (-o)"
    echo
}

## Option defaults
R1=""
outdir=""
repo_dir=""
deinterleave=true

## Get command-line options
while getopts ':i:o:r:dh' flag; do
    case "${flag}" in
    i) R1="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) deinterleave=false ;;
    r) repo_dir="$OPTARG" ;;
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

## Infer the sampleID and define the full output dir
sampleID=$(basename "$R1" | sed 's/_R1.*//')
outdir_full="$outdir"/"$sampleID"

out_mapped="$outdir_full"/mapped_tmp/"$sampleID"
out_unmapped="$outdir_full"/unmapped_tmp/"$sampleID"

R1_mapped="$outdir"/mapped/"$sampleID"_R1_001.fastq.gz
R2_mapped="$outdir"/mapped/"$sampleID"_R2_001.fastq.gz
R1_unmapped="$outdir"/unmapped/"$sampleID"_R1_001.fastq.gz
R2_unmapped="$outdir"/unmapped/"$sampleID"_R2_001.fastq.gz

## Reference FASTA files (to be downloaded)
[[ $repo_dir = "" ]] && repo_dir="$outdir"/"$sampleID"/sortmerna_repo
ref_18s="$repo_dir"/data/rRNA_databases/silva-euk-18s-id95.fasta
ref_28s="$repo_dir"/data/rRNA_databases/silva-euk-28s-id98.fasta

## Report
echo "## R1 FASTQ file:              $R1"
echo "## R2 FASTQ file:              $R2"
echo "## Output dir:                 $outdir_full"
echo "## SortMeRNA repo dir:         $repo_dir"
echo "## 18S reference file:         $ref_18s"
echo "## 28S reference file:         $ref_28s"
echo "## Deinterleave FASTQ files:   $deinterleave"  
echo -e "---------------------------\n"

## Check input
[[ ! -f "$R1" ]] && echo "## ERROR: R1 FASTQ file $R1 not found" >&2 && exit 1
[[ ! -f "$R2" ]] && echo "## ERROR: R2 FASTQ file $R2 not found" >&2 && exit 1

## Make output dirs if needed
mkdir -p "$outdir"/mapped "$outdir"/unmapped "$outdir_full"/mapped_tmp "$outdir_full"/unmapped_tmp


# GET DATABASE FILES -----------------------------------------------------------
## Clone sortmerna repo to get db FASTA files
if [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]]; then
    n_seconds=$(( RANDOM % 50 + 1 ))
    sleep "$n_seconds"s # Sleep for a while so git doesn't error when running this multiple times in parallel
    
    mkdir -p "$repo_dir"
    echo "## Cloning sortmerna repo..."
    [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]] && git clone https://github.com/biocore/sortmerna "$repo_dir"
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
    --aligned "$out_mapped" \
    --other "$out_unmapped" \
    --workdir "$outdir_full" \
    --paired_in \
    --threads "$SLURM_CPUS_PER_TASK"

#?--paired_in Flags the paired-end reads as Aligned, when either of them is Aligned.


# CONVERTING INTERLEAVED FASTQ BACK TO SEPARATED -------------------------------
if [[ "$deinterleave" = true ]]; then
    set +u
    conda deactivate
    source activate "$BBTOOLS_ENV"
    set -u

    echo -e "\n## Deinterleaving R1..."
    reformat.sh \
        in="$out_mapped".fq.gz \
        out1="$R1_mapped" \
        out2="$R2_mapped"

    echo -e "\n## Deinterleaving R2..."
    reformat.sh \
        in="$out_unmapped".fq.gz \
        out1="$R1_unmapped" \
        out2="$R2_unmapped"
    
    echo
else
    mv -v "$out_mapped".fq.gz "$outdir"/mapped
    mv -v "$out_unmapped".fq.gz "$outdir"/unmapped
fi


# HOUSEKEEPING -----------------------------------------------------------------
## Move log files to main dir
mv "$outdir_full"/mapped_tmp/"$sampleID"*log "$outdir"

## Remove temporary files
rm -rv "$outdir_full"/mapped_tmp "$outdir_full"/unmapped_tmp


# QUANTIFY MAPPING SUCCESS -----------------------------------------------------
n_mapped=$(zcat "$R1_mapped" | awk '{ s++ } END{ print s/4 }')
n_unmapped=$(zcat "$R1_unmapped" | awk '{ s++ } END{ print s/4 }')
pct=$(python3 -c "print(round($n_mapped / ($n_unmapped + $n_mapped) * 100, 2))")
echo -e "\nNumber of reads mapped/unmapped, and % mapped:\t$sampleID\t$n_mapped\t$n_unmapped\t$pct"

# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
[[ "$deinterleave" = true ]] && ls -lh "$R1_mapped" "$R2_mapped" "$R1_unmapped" "$R2_unmapped"
echo -e "\n## Done with script sortmerna.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
