#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-bwa-%j.out

## Report
echo -e "\n## Starting script bwa.sh"
date
echo

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bwa-0.7.17

## Bash strict settings
set -euo pipefail

## Command-line args
fq_R1=$1
ref=$2
outdir=$3
readgroup_string=${4-NA}

#fq_R1=data/fastq/AeJap_R1.fastq.gz

## Process args
fq_R2=${fq_R1/_R1/_R2}
sampleID=$(basename "$fq_R1" | sed 's/_R1.*//')  
bam_out=$outdir/"$sampleID".bam
readgroup_arg=""
[[ "$readgroup_string" != NA ]] && readgroup_arg="-R $readgroup_string"

## Check input
[[ ! -f "$fq_R1" ]] && echo "## ERROR: Input R1 FASTQ file $fq_R1 does not exist" >&2 && exit 1
[[ ! -f "$fq_R2" ]] && echo "## ERROR: Input R2 FASTQ file $fq_R2 does not exist" >&2 && exit 1
[[ ! -f "$ref" ]] && echo "## ERROR: Reference FASTA file $ref does not exist" >&2 && exit 1

## Other variables
n_cores=$SLURM_CPUS_PER_TASK

## Make output dir
mkdir -p "$outdir"

## Report
echo "## R1 FASTQ file:           $fq_R1"
echo "## R2 FASTQ file:           $fq_R2"
echo "## Reference FASTA file:    $ref"
echo "## Output BAM file:         $bam_out"
[[ "$readgroup_string" != NA ]] && echo "## Readgroup string:        $readgroup_string"
echo
echo "## Number of cores:         $n_cores"
echo "## SLURM job ID:            $SLURM_JOB_ID"
echo


# MAP WITH BWA -----------------------------------------------------------------
echo -e "## Mapping with bwa mem..."
bwa mem \
    -t "$n_cores" ${readgroup_arg} \
    "$ref" \
    "$fq_R1" "$fq_R1" |
    samtools view -b -h > "$bam_out"

## Report
echo -e "\n## Output BAM file:"
ls -lh "$bam_out"
echo -e "\n## Done with script bwa.sh"
date
echo

################################################################################
## bwa flags:
# -t nr of threads
# -a alignments for single-end / unpaired reads are also output, as secondary alignments
# -M shorter split reads are output as secondary alignments, for Picard compatibility
# -R "@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1"