#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --job-name=STAR-align
#SBATCH --output=slurm-STAR-align-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Align sequences from a FASTQ file to a reference genome with STAR."
  echo
  echo "Syntax: $0 -i <R1-FASTQ-infile> -o <BAM-outdir> -r <ref-index-dir> ..."
  echo
  echo "Required options:"
  echo "   -i FILE      Input gzipped R1 FASTQ file name (The name of the R2 file will be inferred by the script)"
  echo "   -r DIR       Input STAR reference genome index dir (First create index with 'mcic-scripts/rnaseq/star_index.sh')"
  echo "   -o DIR       BAM output dir"
  echo
  echo "Other options:"
  echo "   -a FILE       Reference annotation (GFF/GFF3/GTF) file     [default: no annotation file, but this is not recommended]"
  echo "   -u            Output unmapped reads back as FASTQ file     [default: don't output]"
  echo "   -x            FASTQ files are single-end, not paired-end   [default: FASTQ files are paired-end]"
  echo "   -s            Don't sort the output BAM file               [default: position-sort the BAM file]"
  echo "   -S            Use samtools to sort the output BAM file     [default: Use STAR to sort the BAM file]"
  echo "   -c            Count reads per gene                         [default: don't perform counting]"
  echo "                 NOTE: When using the -c option, you should use a GTF and not a GFF/GFF3"
  echo "   -m INTEGER    Max. nr. of locations a read can map to      [default: 10]"
  echo "   -t INTEGER    Min. intron size                             [default: 21]"
  echo "   -T INTEGER    Max. intron size                             [default: 0 => auto-determined by STAR]"
  echo "   -A STRING     Additional arguments to pass to STAR"
  echo "   -h            Print this help message and exit"
  echo
  echo "Example: $0 -i data/fastq/S01_L001_R1.fastq.gz -o results/star -r refdata/star_index"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
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
count=false
sort=true
samtools_sort=false
output_unmapped=false
paired_end=true
more_args=""

## Parse command-line options
while getopts ':i:o:r:m:t:T:a:A:xuSsch' flag; do
  case "${flag}" in
  i) R1_in="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  o) bam_dir="$OPTARG" ;;
  r) index_dir="$OPTARG" ;;
  m) max_map="$OPTARG" ;;
  t) intron_min="$OPTARG" ;;
  T) intron_max="$OPTARG" ;;
  A) more_args="$OPTARG" ;;
  c) count=true ;;
  s) sort=false ;;
  S) samtools_sort=true ;;
  u) output_unmapped=true ;;
  x) paired_end=false ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SETUP ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env
SAMTOOLS_ENV=/users/PAS0471/jelmer/miniconda3/envs/samtools-env

## Bash strict mode
set -euo pipefail

## Check inputs I
[[ "$R1_in" = "" ]] && echo "## ERROR: Please specify an (R1) input FASTQ file with -i" >&2 && exit 1
[[ "$bam_dir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ "$index_dir" = "" ]] && echo "## ERROR: Please specify a dir with a STAR reference genome index with -r" >&2 && exit 1
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file R1_in (-i) $R1_in does not exist" >&2 && exit 1
[[ ! -d "$index_dir" ]] && echo "## ERROR: Input ref genome dir (-r) $index_dir does not exist" >&2 && exit 1
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && echo "## ERROR: Input annotation file (-a) $gff does not exist" >&2 && exit 1
[[ "$gff" = "" ]] && [[ "$count" = true ]] && echo "## ERROR: Can't do gene counting (-c) without an annotation file (-a)" >&2 && exit 1

## Determine R2 file, output prefix, etc
R1_basename=$(basename "$R1_in" | sed -E 's/.fa?s?t?q.gz//')
if [ "$paired_end" = true ]; then
    ## Determine name of R2 file
    R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}

    sampleID=${R1_basename/"$R1_suffix"/}

    [[ ! -f "$R2_in" ]] && echo "## ERROR: Input file R2_in $R2_in does not exist" >&2 && exit 1
    [[ "$R1_in" = "$R2_in" ]] && echo "## ERROR: Input file R1 is the same as R2" >&2 && exit 1
else
    sampleID="$R1_basename"
    R2_in=""
fi
outfile_prefix="$bam_dir/$sampleID"_

## Other output dirs
starlog_dir="$bam_dir"/star_logs

## If a GFF/GTF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then

    if [[ "$gff" =~ .*\.gff3? ]]; then
        annot_format=gff
        GFF_FORMAT="--sjdbGTFtagExonParentTranscript Parent"
    elif [[ "$gff" =~ .*\.gff3? ]]; then
        annot_format=gtf
        GFF_FORMAT="--sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFtagExonParentGene gene_id"
    else
        echo "## ERROR: Unknown annotation file format" && exit 1
    fi

    if [ "$count" = true ]; then
        if [ "$annot_format" = "gff" ]; then
            # Better to use GTF for cournting https://groups.google.com/g/rna-star/c/M0q8M5FscA4
            echo "## ERROR: Please convert your GFF to a GTF, use the script 'mcic-scripts/convert/gff2gtf.sh'"
        fi
        # add `--quantMode GeneCounts` so as to do gene counting 
        annot_arg="--sjdbGTFfile $gff $GFF_FORMAT --quantMode GeneCounts"
    else
        annot_arg="--sjdbGTFfile $gff $GFF_FORMAT"
    fi
else
    annot_arg=""
fi

## Sorted output or not
if [[ "$sort" = true && "$samtools_sort" = false ]]; then
    output_arg="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 100"
else
    output_arg="--outSAMtype BAM Unsorted"
fi

## Output unmapped reads in a FASTQ file
if [ "$output_unmapped" = true ]; then
    unmapped_arg="--outReadsUnmapped Fastx"
    unmapped_dir="$bam_dir"/unmapped
    mkdir -p "$unmapped_dir"
else
    unmapped_arg=""
fi

## Report
echo
echo "## Starting script star_align.sh"
date
echo
echo "## Input (R1) FASTQ file:                        $R1_in"
echo "## Input STAR genome index dir:                  $index_dir"
[[ "$gff" != "" ]] && echo "## Input GFF file:                               $gff"
echo "## Output BAM dir:                               $bam_dir"
echo
echo "## Output unmapped reads as FASTQ:               $output_unmapped"
echo "## Also perform read counting:                   $count"
echo "## Max nr of alignments for a read:              $max_map"  # If this nr is exceeded, read is considered unmapped
echo "## Minimum intron size:                          $intron_min"
echo "## Maximum intron size (0 => STAR default):      $intron_max"
echo "## Sort the BAM file:                            $sort"
echo "## Sort the BAM file with samtools:              $samtools_sort"
[[ "$more_args" != "" ]] && echo "## Additional args to pass to STAR:              $more_args"
echo
echo "## Sample ID (as inferred by the script):        $sampleID"
[[ "$paired_end" != "" ]] && echo "## R2 FASTQ file (as inferred by the script):    $R2_in"
[[ "$gff" != "" ]] && echo "## Annotation arg for STAR:                      $annot_arg"
echo "## Output arg for STAR:                          $output_arg"
echo -e "------------------------\n"

## Create output dirs if needed
mkdir -p "$bam_dir" "$starlog_dir"


# ALIGN ------------------------------------------------------------------------
echo "## Aligning reads with STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
     --genomeDir "$index_dir" \
     --readFilesIn "$R1_in" "$R2_in" \
     --readFilesCommand zcat \
     --outFilterMultimapNmax $max_map \
     --alignIntronMin $intron_min --alignIntronMax $intron_max \
     --outFileNamePrefix "$outfile_prefix" \
     $unmapped_arg $output_arg $annot_arg $more_args


# SORT WITH SAMTOOLS SORT ------------------------------------------------------
## STAR may fail to sort on some large BAM files, or BAM files with a lot of
## reads mapping to similar positions. In that case, could use samtools sort. 

if [[ "$sort" = true && "$samtools_sort" = "true" ]]; then

    echo -e "\n## Sorting the BAM file with samtools sort..."
    
    bam_unsorted="$outfile_prefix"Aligned.out.bam
    bam_sorted="$outfile_prefix"Aligned.sortedByCoord.out.bam
    
    source activate "$SAMTOOLS_ENV"
    samtools sort "$bam_unsorted" > "$bam_sorted"
fi


# ORGANIZE STAR OUTPUT ---------------------------------------------------------
if [ "$output_unmapped" = true ]; then
    echo -e "\n## Moving, renaming and zipping unmapped FASTQ files...."
    for oldpath in "$outfile_prefix"*Unmapped.out.mate*; do
        oldname=$(basename "$oldpath")
        newname=$(echo "$oldname" | sed -E s'/_Unmapped.out.mate([12])/_R\1.fastq.gz/')
        newpath="$unmapped_dir"/"$newname"

        echo "## Moving $oldpath to $newpath..."

        #> The unmapped FASTQ files output by STAR have a weird format with "0:N" for R1 reads (instead of "1:N")
        #> and "1:N" (instead of "2:N") for R2 reads, which Trinity doesn't accept. The code below will fix that:
        [[ "$newpath" = *R1.fastq.gz ]] && sed -E 's/(^@.*) 0:N: (.*)/\1 1:N: \2/' "$oldpath" | gzip -f > "$newpath"
        [[ "$newpath" = *R2.fastq.gz ]] && sed -E 's/(^@.*) 1:N: (.*)/\1 2:N: \2/' "$oldpath" | gzip -f > "$newpath"

        ## Remove old file
        rm -v "$oldpath"
    done
fi

## Move STAR log files
echo -e "\n## Moving STAR log files...."
mv -v "$bam_dir"/"$sampleID"*out "$starlog_dir"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "## Listing output BAM file(s):"
ls -lh "$outfile_prefix"*bam

if [ "$output_unmapped" = true ]; then
    echo "## Listing unmapped FASTQ files:"
    ls -lh "$unmapped_dir/$sampleID"*fastq.gz
fi

if [ "$count" = true ]; then
    echo -e "\n## Listing output gene count file:"
    ls -lh "$bam_dir"/"$sampleID"*ReadsPerGene.out.tab
fi

echo -e "\n## Done with script star_align.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
