#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --job-name=STAR-align
#SBATCH --output=slurm-STAR-align-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                             $0"
    echo "                 ALIGN RNASEQ READS TO A GENOME WITH STAR"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "sbatch $0 -i <R1-FASTQ> -r <ref-index-dir> -o <outdir> [...]"
    echo
    echo "REQUIRED OPTIONS:"
    echo "   -i FILE       Input gzipped R1 FASTQ file path (The name of the R2 file will be inferred by the script)"
    echo "   -r DIR        Input STAR reference genome index dir (First create index with 'mcic-scripts/rnaseq/star_index.sh')"
    echo "   -o DIR        BAM output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "   -a FILE       Reference annotation (GFF/GFF3/GTF) file (GTF preferred)  [default: no annotation file, but this is not recommended]"
    echo "   -u            Output unmapped reads back as FASTQ file             [default: don't output]"
    echo "   -P            FASTQ files are single-end, not paired-end           [default: paired-end]"
    echo "   -s            Don't sort the output BAM file                       [default: position-sort the BAM file]"
    echo "   -S            Use samtools to sort the output BAM file             [default: Use STAR to sort the BAM file]"
    echo "   -c            Count reads per gene                                 [default: don't perform counting]"
    echo "                 NOTE: When using the -c option, you should use a GTF and not a GFF/GFF3"
    echo "   -m INTEGER    Max. nr. of locations a read can map to              [default: 10]"
    echo "   -t INTEGER    Min. intron size                                     [default: 21 (also the STAR default)]"
    echo "   -T INTEGER    Max. intron size                                     [default: 0 => auto-determined by STAR]"
    echo "   -A STRING     Quoted string with additional arguments to pass to STAR"
    echo
    echo
    echo "UTILITY OPTIONS:"
    echo "   -h            Print this help message and exit"
    echo "   -x            Turn on debugging mode: print all code in the script"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "sbatch $0 -i data/fastq/S01_R1.fastq.gz -o results/star -r refdata/star_index"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/star-2.7.10a  # NOTE: This env includes samtools
}

## Print version
Print_version() {
    Load_software
    STAR --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
max_map=10
intron_min=21             # STAR default, too
intron_max=0              # => auto-determined; STAR default, too
count=false
sort=true
samtools_sort=false
output_unmapped=false && unmapped_arg=""
paired_end=true
debug=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
R1_in=""
outdir=""
index_dir=""
annot="" && annot_arg=""
more_args=""

## Parse command-line options
while getopts ':i:o:r:m:t:T:a:A:PxuSsch' flag; do
    case "${flag}" in
        i) R1_in="$OPTARG" ;;
        a) annot="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        r) index_dir="$OPTARG" ;;
        m) max_map="$OPTARG" ;;
        t) intron_min="$OPTARG" ;;
        T) intron_max="$OPTARG" ;;
        A) more_args="$OPTARG" ;;
        c) count=true ;;
        s) sort=false ;;
        S) samtools_sort=true ;;
        u) output_unmapped=true ;;
        P) paired_end=false ;;
        x) debug=true ;;
        h) Print_help && exit 0 ;;
        \?) Die "Invalid option" ;;
        :) Die "Option -$OPTARG requires an argument" ;;
    esac
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
Load_software

## Bash strict mode
set -euo pipefail

## Check inputs I
[[ "$R1_in" = "" ]] && Die "Please specify an (R1) input FASTQ file with -i"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ "$index_dir" = "" ]] && Die "Please specify a dir with a STAR reference genome index with -r"
[[ ! -f "$R1_in" ]] && Die "Input file R1_in (-i) $R1_in does not exist"
[[ ! -d "$index_dir" ]] && Die "Input ref genome dir (-r) $index_dir does not exist"
[[ "$annot" != "" ]] && [[ ! -f "$annot" ]] && Die "Input annotation file (-a) $annot does not exist"
[[ "$annot" = "" ]] && [[ "$count" = true ]] && Die "Need an annotation file (-a) for gene counting (-c)"

## Determine R2 file, output prefix, etc
R1_basename=$(basename "$R1_in" | sed -E 's/.fa?s?t?q.gz//')

if [ "$paired_end" = true ]; then
    
    ## Determine name of R2 file
    R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}

    sampleID=${R1_basename/"$R1_suffix"/}

    [[ ! -f "$R2_in" ]] && Die "Input file R2_in $R2_in does not exist"
    [[ "$R1_in" = "$R2_in" ]] && Die "Input file R1 is the same as R2"
else
    
    sampleID="$R1_basename"
    R2_in=""
fi

outfile_prefix="$outdir/$sampleID"_

## Other output dirs
starlog_dir="$outdir"/star_logs

## If a GFF/GTF file is provided, build the appropriate argument for STAR
if [ "$annot" != "" ]; then

    if [[ "$annot" =~ .*\.gff3? ]]; then
        annot_format=gff
        annot_tags="--sjdbGTFtagExonParentTranscript Parent"
    elif [[ "$annot" =~ .*\.gtf ]]; then
        annot_format=gtf
        annot_tags="--sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFtagExonParentGene gene_id"
    else
        Die "Unknown annotation file format"
    fi

    if [ "$count" = true ]; then
        if [ "$annot_format" = "gff" ]; then
            # Better to use GTF for counting https://groups.google.com/g/rna-star/c/M0q8M5FscA4
            Die "Please convert your GFF to a GTF, use the script 'mcic-scripts/convert/gff2gtf.sh'"
        fi
        # add `--quantMode GeneCounts` so as to do gene counting 
        annot_arg="--sjdbGTFfile $annot $annot_tags --quantMode GeneCounts"
    else
        annot_arg="--sjdbGTFfile $annot $annot_tags"
    fi

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
    unmapped_dir="$outdir"/unmapped
    mkdir -p "$unmapped_dir"
fi

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT STAR_ALIGN.SH"
date
echo "=========================================================================="
echo "Input (R1) FASTQ file:                        $R1_in"
echo "Input STAR genome index dir:                  $index_dir"
[[ "$annot" != "" ]] && echo "Input annotation file:                        $annot"
echo "Output BAM dir:                               $outdir"
echo
echo "Output unmapped reads as FASTQ:               $output_unmapped"
echo "Also perform read counting:                   $count"
echo "Max nr of alignments for a read:              $max_map"  # If this nr is exceeded, read is considered unmapped
echo "Minimum intron size:                          $intron_min"
echo "Maximum intron size (0 => STAR default):      $intron_max"
echo "Sort the BAM file:                            $sort"
echo "Sort the BAM file with samtools:              $samtools_sort"
echo "Are FASTQ reads paired-end?                   $paired_end"
[[ "$more_args" != "" ]] && echo "Additional args to pass to STAR:              $more_args"
echo
echo "Sample ID (as inferred by the script):        $sampleID"
[[ "$paired_end" = true ]] && echo "R2 FASTQ file (as inferred by the script):    $R2_in"
[[ "$annot" != "" ]] && echo "Annotation arg for STAR:                      $annot_arg"
echo "Output arg for STAR:                          $output_arg"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
mkdir -p "$outdir"/logs "$starlog_dir"

## Run STAR
echo "## Aligning reads with STAR...."
STAR --runThreadN "$SLURM_CPUS_ON_NODE" \
    --genomeDir "$index_dir" \
    --readFilesIn "$R1_in" "$R2_in" \
    --readFilesCommand zcat \
    --outFilterMultimapNmax $max_map \
    --alignIntronMin $intron_min \
    --alignIntronMax $intron_max \
    --outFileNamePrefix "$outfile_prefix" \
    $unmapped_arg $output_arg $annot_arg $more_args


## Sort with samtools sort
## STAR may fail to sort on some large BAM files, or BAM files with a lot of
## reads mapping to similar positions. In that case, could use samtools sort. 
if [[ "$sort" = true && "$samtools_sort" = "true" ]]; then

    echo -e "\n## Sorting the BAM file with samtools sort..."
    
    bam_unsorted="$outfile_prefix"Aligned.out.bam
    bam_sorted="$outfile_prefix"Aligned.sortedByCoord.out.bam
    
    samtools sort "$bam_unsorted" > "$bam_sorted"
fi

## Organize STAR output
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
mv -v "$outdir"/"$sampleID"*out "$starlog_dir"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo "========================================================================="
echo -e "## Listing output BAM file(s):"
ls -lh "$outfile_prefix"*bam

if [ "$output_unmapped" = true ]; then
    echo "## Listing unmapped FASTQ files:"
    ls -lh "$unmapped_dir/$sampleID"*fastq.gz
fi

if [ "$count" = true ]; then
    echo -e "\n## Listing output gene count file:"
    ls -lh "$outdir"/"$sampleID"*ReadsPerGene.out.tab
fi

echo
echo "## Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
echo "## Done with script"
date
