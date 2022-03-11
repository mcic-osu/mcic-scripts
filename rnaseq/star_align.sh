#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --job-name=STAR_align
#SBATCH --output=slurm-STAR-align-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help
Help() {
  echo
  echo "## $0: Align sequences from a FASTQ file to a reference genome with STAR."
  echo
  echo "## Syntax: $0 -i <R1-FASTQ-infile> -o <BAM-outdir> -r <ref-index-dir> ..."
  echo
  echo "## Required options:"
  echo "##    -i STR    R1 FASTQ input file name (The name of the R2 file will be inferred by the script.)"
  echo "##    -r STR    STAR reference genome index dir"
  echo "##    -o STR    BAM output dir"
  echo
  echo "## Other options:"
  echo "##    -a STR    Reference annotation (GFF/GTF) file (default: no GFF/GTF, but this is not recommended)"
  echo "##    -c        Count reads per gene (default: don't count)"
  echo "##    -m INT    Max. number of locations a read can map to, before being considered unmapped (default: 10)"
  echo "##    -t INT    Min. intron size (default: 21)"
  echo "##    -T INT    Max. intron size (default: 0 => determined by STAR)"
  echo "##    -A STR    Additional arguments to pass to STAR"
  echo "##    -h        Print this help message"
  echo
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
do_count=false

## Parse command-line options
while getopts ':i:o:r:m:t:T:a:A:ch' flag; do
  case "${flag}" in
  i) R1_in="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  o) bam_dir="$OPTARG" ;;
  r) index_dir="$OPTARG" ;;
  m) max_map="$OPTARG" ;;
  t) intron_min="$OPTARG" ;;
  T) intron_max="$OPTARG" ;;
  A) more_args="$OPTARG" ;;
  c) do_count=true ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SETUP ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env

## Bash strict mode
set -euo pipefail

## Hardcoded parameters
GFF_FORMAT="--sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent"

## Process parameters
R2_in=${R1_in/_R1_/_R2_}                            # R2 FASTQ file
sample_id=$(basename "$R1_in" | sed 's/_R1_.*//')   # Sample ID - for output files
unmapped_dir="$bam_dir"/unmapped                    # STAR output other than bam files
starlog_dir="$bam_dir"/star_logs                    # STAR output other than bam files

## If a GFF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then
    if [ "$do_count" = true ]; then
        # add `--quantMode GeneCounts` so as to do gene counting 
        gff_arg="--sjdbGTFfile $gff $GFF_FORMAT --quantMode GeneCounts"
    else
        gff_arg="--sjdbGTFfile $gff $GFF_FORMAT"
    fi
else
    gff_arg=""
fi

## Report
echo "## Starting script star_align.sh"
date
echo

## Check inputs
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file R1_in (-i) $R1_in does not exist" >&2 && exit 1
[[ ! -f "$R2_in" ]] && echo "## ERROR: Input file R2_in $R2_in does not exist" >&2 && exit 1
[[ ! -d "$index_dir" ]] && echo "## ERROR: Input ref genome dir (-r) $index_dir does not exist" >&2 && exit 1
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && echo "## ERROR: Input annotation file (-a) $gff does not exist" >&2 && exit 1
[[ "$gff" = "" ]] && [[ "$do_count" = true ]] && echo "## ERROR: Can't do gene counting (-c) without an annotation file (-a)" >&2 && exit 1

## Report
echo "## Input R1 FASTQ file:                          $R1_in"
echo "## Input STAR genome index dir:                  $index_dir"
[[ "$gff" != "" ]] && echo "## Input GFF file:                               $gff"
echo "## Output BAM dir:                               $bam_dir"
echo
echo "## Also perform read counting:                   $do_count"
echo "## Max nr of alignments for a read:              $max_map"  # If this nr is exceeded, read is considered unmapped
echo "## Min intron size:                              $intron_min"
echo "## Max intron size (0 => STAR default):          $intron_max"
echo "## Additional args to pass to STAR:              $more_args"
echo
echo "## Sample ID (as inferred by the script):        $sample_id"
echo "## R2 FASTQ file (as inferred by the script):    $R2_in"
echo -e "------------------------\n"

## Create output dirs if needed
mkdir -p "$bam_dir"
mkdir -p "$starlog_dir"
mkdir -p "$unmapped_dir"


# ALIGN ------------------------------------------------------------------------
echo "## Aligning reads with STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
   --genomeDir "$index_dir" \
   --readFilesIn "$R1_in" "$R2_in" \
   --readFilesCommand zcat \
   --outFilterMultimapNmax $max_map \
   --alignIntronMin $intron_min --alignIntronMax $intron_max \
   --outFileNamePrefix "$bam_dir/$sample_id"_ \
   --outSAMtype BAM SortedByCoordinate \
   --outBAMsortingBinsN 100 \
   --outReadsUnmapped Fastx $gff_arg $more_args


# ORGANIZE STAR OUTPUT ---------------------------------------------------------
## Move files with unmapped reads
echo -e "\n## Moving, renaming and zipping unmapped FASTQ files...."
for oldpath in "$bam_dir/$sample_id"_*Unmapped.out.mate*; do
    oldname=$(basename "$oldpath")
    newname=$(echo "$oldname" | sed -E s'/_Unmapped.out.mate([12])/_R\1.fastq.gz/')
    newpath="$unmapped_dir"/"$newname"

    #> The unmapped FASTQ files output by STAR have a weird format with "0:N" for R1 reads (instead of "1:N")
    #> and "1:N" (instead of "2:N") for R2 reads, which Trinity doesn't accept. The code below will fix that:
    [[ "$newpath" = *R1.fastq.gz ]] && sed -E 's/(^@.*) 0:N: (.*)/\1 1:N: \2/' "$oldpath" | gzip -f > "$newpath"
    [[ "$newpath" = *R2.fastq.gz ]] && sed -E 's/(^@.*) 1:N: (.*)/\1 2:N: \2/' "$oldpath" | gzip -f > "$newpath"
done

## Move STAR log files
echo -e "\n## Moving STAR log files...."
mv -v "$bam_dir"/"$sample_id"*out "$starlog_dir"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"

echo "## Listing unmapped FASTQ files:"
ls -lh "$unmapped_dir/$sample_id"*fastq.gz

echo -e "\n## Listing output BAM file:"
ls -lh "$bam_dir/$sample_id"_*bam

if [ "$do_count" = true ]; then
    echo -e "\n## Listing output gene count file:"
    ls -lh "$bam_dir"/"$sample_id"*ReadsPerGene.out.tab
fi

echo -e "\n## Done with script star_align.sh"
date
