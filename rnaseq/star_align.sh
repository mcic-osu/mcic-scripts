#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --job-name=STAR_align
#SBATCH --output=slurm-STAR-align-%j.out


# SETUP ---------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env
conda activate --stack /users/PAS0471/jelmer/.conda/envs/samtools-env

## Bash strict mode
set -euo pipefail

## Help
Help() {
  echo
  echo "## $0: Align sequences from a FASTQ file to a reference genome with STAR."
  echo
  echo "## Syntax: $0 -i <R1 FASTQ input file> -o <BAM output dir> -r <ref genome index dir> [ -m <max multi-map> ] [-h]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   R1 FASTQ input file (REQUIRED; note that the name of the R2 file will be inferred by the script.)"
  echo "## -o STR   BAM output dir (REQUIRED)"
  echo "## -r STR   Reference index dir (REQUIRED)"
  echo "## -g STR   Reference GFF file (REQUIRED)"
  echo "## -m INT   Max. number of locations a read can map to, before being considered unmapped (default: 10)"
  echo "## -t INT   Min. intron size (default: 21)"
  echo "## -T INT   Max. intron size (default: 0 => determined by STAR)"
  echo "## Example: $0 -i data/fastq/S01_L001_R1.fastq.gz -o results/star -r refdata/star_index"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
R1_in=""
bam_dir=""
index_dir=""
gff=""
max_map=10
intron_min=21
intron_max=0

## Parse command-line options
while getopts ':i:o:r:m:t:T:g:h' flag; do
  case "${flag}" in
  i) R1_in="$OPTARG" ;;
  g) gff="$OPTARG" ;;
  o) bam_dir="$OPTARG" ;;
  r) index_dir="$OPTARG" ;;
  m) max_map="$OPTARG" ;;
  t) intron_min="$OPTARG" ;;
  T) intron_max="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Process options
R2_in=${R1_in/_R1_/_R2_}
sample_id=$(basename "$R1_in" | sed 's/_R1_.*//')

flagstat_dir="$bam_dir"/flagstat
flagstat_out="$flagstat_dir"/"$sample_id"_flagstat.txt

unmapped_dir="$bam_dir"/unmapped
starlog_dir="$bam_dir"/star_logs

## Report
echo "## Starting script star_align.sh"
date
echo "## Input R1 FASTQ file:                 $R1_in"
echo "## Input STAR genome index dir:         $index_dir"
echo "## Input GFF file:                      $gff"
echo "## Output BAM dir:                      $bam_dir"
echo
echo "## Max nr of alignments for a read:     $max_map"  # If this nr is exceeded, read is considered unmapped
echo "## Min intron size:                     $intron_min"
echo "## Max intron size (0 => STAR default): $intron_max"
echo
echo "## Sample ID (inferred):                $sample_id"
echo "## R2 FASTQ file (input - inferred):    $R2_in"
echo -e "------------------------\n"

## Check inputs
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file $R1_in does not exist" >&2 && exit 1
[[ ! -f "$R2_in" ]] && echo "## ERROR: Input file $R2_in does not exist" >&2 && exit 1
[[ ! -f "$gff" ]] && echo "## ERROR: Input file $gff does not exist" >&2 && exit 1
[[ ! -d "$index_dir" ]] && echo "## ERROR: Input dir $index_dir does not exist" >&2 && exit 1

## Create output dirs if needed
mkdir -p "$bam_dir"
mkdir -p "$flagstat_dir"
mkdir -p "$starlog_dir"
mkdir -p "$unmapped_dir"


# ALIGN ------------------------------------------------------------------------
echo "## Aligning reads with STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
   --genomeDir "$index_dir" \
   --sjdbGTFfile "$gff" \
   --readFilesIn "$R1_in" "$R2_in" \
   --readFilesCommand zcat \
   --outFilterMultimapNmax $max_map \
   --alignIntronMin $intron_min --alignIntronMax $intron_max \
   --outFileNamePrefix "$bam_dir/$sample_id"_ \
   --outSAMtype BAM SortedByCoordinate \
   --outBAMsortingBinsN 100 \
   --outReadsUnmapped Fastx \
   --quantMode GeneCounts


# ORGANIZE STAR OUTPUT ---------------------------------------------------------
## Move files with unmapped reads
echo -e "\n## Moving, renaming and zipping unmapped FASTQ files...."
for oldpath in "$bam_dir/$sample_id"_*Unmapped.out.mate*; do
    oldname=$(basename "$oldpath")
    newname=$(echo "$oldname" | sed -E s'/_Unmapped.out.mate([12])/_R\1.fastq/')
    newpath="$unmapped_dir"/"$newname"

    echo "## Old filename: $oldpath"
    echo "## New filename: $newpath"

    mv -v "$oldpath" "$newpath"
    gzip -f "$newpath"

    echo
done

## Move STAR log files
echo -e "\n## Moving STAR log files...."
mv -v "$bam_dir"/"$sample_id"*out "$bam_dir"/"$sample_id"*tab "$starlog_dir"


# RUN SAMTOOLS FLAGSTAT --------------------------------------------------------
echo -e "\n## Running samtools flagstat...."
samtools flagstat "$bam_dir/$sample_id"_*bam > "$flagstat_out"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing output BAM file:"
ls -lh "$bam_dir/$sample_id"_*bam
echo "## Listing unmapped FASTQ files:"
ls -lh "$unmapped_dir/$sample_id"*fastq.gz
echo "## Listing flagstat file:"
ls -lh "$flagstat_out"

echo -e "\n## Done with script star_align.sh"
date
