#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-bwa-%j.out


# SETUP ------------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Map FASTQ reads to a reference genome with BWA."
  echo
  echo "## Syntax: $0 -i <R1-FASTQ-file> -r <ref-genome> -o <output-dir> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING        Input (R1) FASTQ file (if paired-end, name of R2 file will be inferred)"
  echo "## -o STRING        Output directory"
  echo "## -r STRING        Reference genome FASTA file (with BWA indices in same dir)"
  echo
  echo "## Other options:"
  echo "## -s               Sequences are single-end -- don't infer name of R2 file [default: paired-end sequences]"
  echo "## -g STRING        Readgroup string to be added"
  echo "## -h               Print this help message and exit"
  echo
  echo "## Example command:"
  echo "## $0 -i data/fastq/A1_R1.fastq.gz -o results/bwa -r refdata/genome.fa"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
fq_R1=""
ref=""
outdir=""
readgroup_string=""
single_end=false

## Get command-line options
while getopts 'i:o:r:g:sh' flag; do
    case "${flag}" in
    i) fq_R1="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    r) ref="$OPTARG" ;;
    g) readgroup_arg="$OPTARG" ;;
    s) single_end=true ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script bwa.sh"
date
echo

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bwa-0.7.17
SAMTOOLS_ENV=/fs/ess/PAS0471/jelmer/conda/samtools # Loaded later

## Bash strict settings
set -euo pipefail

## Process parameters - input file arg
if [[ "$single_end" = false ]]; then
    fq_R2=${fq_R1/_R1/_R2}
    fq_arg="$fq_R1 $fq_R2"
else
    fq_arg="$fq_R1"
fi

## Process parameters -- output files
sampleID=$(basename "$fq_R1" | sed -E 's/_R?1[\._].*//')  
bam_out=$outdir/"$sampleID".bam
flagstat_out=$outdir/"$sampleID".flagstat

## Process parameters -- readgroup arg
readgroup_arg=""
[[ "$readgroup_string" != "" ]] && readgroup_arg="-R $readgroup_string"

## Check input
[[ ! -f "$fq_R1" ]] && echo "## ERROR: Input R1 FASTQ file $fq_R1 does not exist" >&2 && exit 1
[[ "$single_end" = false ]] && [[ ! -f "$fq_R2" ]] && echo "## ERROR: Input R2 FASTQ file $fq_R2 does not exist" >&2 && exit 1
[[ ! -f "$ref" ]] && echo "## ERROR: Reference FASTA file $ref does not exist" >&2 && exit 1

## Other variables
n_cores=$SLURM_CPUS_PER_TASK

## Make output dir
mkdir -p "$outdir"

## Report
echo "## R1 FASTQ file:             $fq_R1"
[[ "$single_end" = false ]] && echo "## R2 FASTQ file:             $fq_R2"
echo "## Reference FASTA file:      $ref"
echo "## Output BAM file:           $bam_out"
[[ "$readgroup_string" != NA ]] && echo "## Readgroup string:        $readgroup_string"
echo
echo "## Number of cores:           $n_cores"
echo "## SLURM job ID:              $SLURM_JOB_ID"
echo


# MAIN -------------------------------------------------------------------------
## Map
echo "## Mapping with bwa mem..."
bwa mem \
    -t "$n_cores" ${readgroup_arg} \
    "$ref" \
    ${fq_arg} |
    samtools view -b -h > "$bam_out"

## Get mapping stats
echo -e "\n## Running samtools flagstat..."
source activate $SAMTOOLS_ENV
samtools flagstat "$bam_out" > "$flagstat_out"


# WRAP UP ----------------------------------------------------------------------
## Report
echo -e "\n## Output files:"
ls -lh "$bam_out" "$flagstat_out"
echo -e "\n## Done with script bwa.sh"
date
echo

################################################################################
## bwa flags:
# -t nr of threads
# -a alignments for single-end / unpaired reads are also output, as secondary alignments
# -M shorter split reads are output as secondary alignments, for Picard compatibility
# -R "@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1"