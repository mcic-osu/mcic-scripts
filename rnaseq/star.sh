#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=star-align
#SBATCH --output=slurm-star_align-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                             $0"
    echo "                 ALIGN RNASEQ READS TO A GENOME WITH STAR"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 [ -i <R1-FASTQ> --fofn <FOFN> ] -r <ref-index-dir> -o <outdir> [...]"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -r/--index_dir      <dir>   Input STAR reference genome index dir (First create index with 'mcic-scripts/rnaseq/star_index.sh')"
    echo "  -o/--outdir         <dir>   BAM output dir"
    echo "To specify the input reads, use one of the following to options:"
    echo "  -i/--R1/--reads     <file>  Input gzipped (R1) FASTQ file path (If R1, name of R2 will be inferred unless using '--single_end')"
    echo "  --fofn              <file>  A File of File Names (FOFN), with one line per input file"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --annot             <file>  Reference annotation (GFF/GFF3/GTF) file (GTF preferred)  [default: no annotation file, but this is not recommended]"
    echo "  --output_unmapped           Output unmapped reads back as FASTQ file    [default: don't output]"
    echo "  --single_end                FASTQ files are single-end, not paired-end  [default: paired-end]"
    echo "  --no_sort                   Don't sort the output BAM file              [default: position-sort the BAM file]"
    echo "  --samtools_sort             Use samtools to sort the output BAM file    [default: Use STAR to sort the BAM file]"
    echo "  --index_bam                 Index the output BAM file with samtools     [default: don't index]"
    echo "  --count                     Count reads per gene                        [default: don't perform counting]"
    echo "                              NOTE: When using this flag, provide a GTF and not a GFF/GFF3 file"
    echo "  --max_map           <int>   Max. nr. of locations a read can map to     [default: 10]"
    echo "  --intron_min        <int>   Min. intron size                            [default: 21 (also the STAR default)]"
    echo "  --intron_max        <int>   Max. intron size                            [default: 0 => auto-determined by STAR]"
    echo "  --more_args         <str>   Quoted string with additional argument(s) to pass to STAR"
    echo
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for START and exit"
    echo "  -v/--version                Print the version of STAR and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/S01_R1.fastq.gz -o results/star -r refdata/star_index"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - https://github.com/alexdobin/STAR"
    echo "  - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/star-2.7.10a  # NOTE: This env includes samtools
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    STAR --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    STAR --help
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
}

# Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Account (project):    $SLURM_JOB_ACCOUNT"
    echo "Job ID:               $SLURM_JOB_ID"
    echo "Job name:             $SLURM_JOB_NAME"
    echo "Memory (per node):    $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

# Set the number of threads/CPUs
Set_threads() {
    set +u
    if [[ "$slurm" = true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            threads="$SLURM_NTASKS"
        else 
            echo "WARNING: Can't detect nr of threads, setting to 1"
            threads=1
        fi
    else
        threads=1
    fi
    set -u
}

# Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo >&2
    echo "=====================================================================" >&2
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option" >&2
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'" >&2
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    echo -e "\nEXITING..." >&2
    echo "=====================================================================" >&2
    echo >&2
    exit 1
}



# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
max_map=10                # If this nr is exceeded, read is considered unmapped
intron_min=21             # STAR default, too
intron_max=0              # => auto-determined; STAR default, too
count=false
sort_bam=true
samtools_sort=false
index_bam=false
output_unmapped=false && unmapped_arg=""
paired_end=true

debug=false
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
R1_in=""
R2_in=""
fofn=""
declare -a infiles
outdir=""
index_dir=""
annot="" && annot_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 | --reads )   shift && R1_in=$1 ;;
        --fofn )                shift && fofn=$1 ;;
        -r | --index_dir )      shift && index_dir=$1 ;;
        -a | --annot )          shift && annot=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --max_map )             shift && max_map=$1 ;;
        --intron_min )          shift && intron_min=$1 ;;
        --intron_max )          shift && intron_max=$1 ;;
        --output_unmapped )     shift && output_unmapped=$1 ;;
        --count )               count=true ;;
        --index_bam )           index_bam=true ;;
        --no_sort )             sort_bam=false ;;
        --samtools_sort )       samtools_sort=true ;;
        --single_end )          paired_end=false ;;
        --more_args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0;;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash strict mode
set -euo pipefail

# Load software
Load_software
Set_threads

# Check inputs I
[[ "$R1_in" = "" && "$fofn" = "" ]] && Die "Please specify input FASTQ file(s) with -i or --fofn"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ "$index_dir" = "" ]] && Die "Please specify a dir with a STAR reference genome index with -r"
[[ ! -d "$index_dir" ]] && Die "Input ref genome dir (-r) $index_dir does not exist"
[[ "$annot" != "" ]] && [[ ! -f "$annot" ]] && Die "Input annotation file (-a) $annot does not exist"
[[ "$annot" = "" ]] && [[ "$count" = true ]] && Die "Need an annotation file (-a) for gene counting (-c)"

# Input files via FOFN
if  [[ "$fofn" != "" ]]; then
    mapfile -t infiles <"$fofn"
    R1_in=${infiles[0]}
    [[ ${#infiles[@]} -eq 2 ]] && R2_in=${infiles[1]}
    [[ ${#infiles[@]} -gt 2 ]] && Die "FOFN should contain 1 or 2 filenames, not ${#infiles[@]}"
fi

# Determine R2 file, output prefix, etc
if [[ "$R1_in" != "" ]]; then
    R1_basename=$(basename "$R1_in" | sed -E 's/.fa?s?t?q.gz//')
    
    if [ "$paired_end" = true ]; then
        R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
        sampleID=${R1_basename/"$R1_suffix"/}

        # Determine name of R2 file
        if [[ "$R2_in" = "" ]]; then
            R2_suffix=${R1_suffix/1/2}
            R2_in=${R1_in/$R1_suffix/$R2_suffix}
        fi
        
        [[ ! -f "$R1_in" ]] && Die "Input file R1_in $R1_in does not exist"
        [[ ! -f "$R2_in" ]] && Die "Input file R2_in $R2_in does not exist"
        [[ "$R1_in" = "$R2_in" ]] && Die "Input file R1 is the same as R2"
    
    else
        sampleID="$R1_basename"
    fi
fi

# Define other outputs
outfile_prefix="$outdir/$sampleID"_
starlog_dir="$outdir"/star_logs
final_bam="$outdir"/bam/"$sampleID".bam

# If a GFF/GTF file is provided, build the appropriate argument for STAR
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

# Sorted output or not
if [[ "$sort_bam" = true && "$samtools_sort" = false ]]; then
    output_arg="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 100"
else
    output_arg="--outSAMtype BAM Unsorted"
fi

# Output unmapped reads in a FASTQ file
if [ "$output_unmapped" = true ]; then
    unmapped_arg="--outReadsUnmapped Fastx"
    unmapped_dir="$outdir"/unmapped
    mkdir -p "$unmapped_dir"
fi

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT STAR_ALIGN.SH"
date
echo "=========================================================================="
echo "Output BAM dir:                               $outdir"
echo "Input R1 FASTQ file:                          $R1_in"
[[ "$R2_in" != "" ]] && echo "Input R2 FASTQ file:  $R2_in"
echo "Are FASTQ reads paired-end?                   $paired_end"
echo "Input STAR genome index dir:                  $index_dir"
[[ "$annot" != "" ]] && echo "Input annotation file:                        $annot"
echo "Output unmapped reads as FASTQ:               $output_unmapped"
echo "Also perform read counting:                   $count"
echo "Max nr of alignments for a read:              $max_map"
echo "Minimum intron size:                          $intron_min"
echo "Maximum intron size (0 => STAR default):      $intron_max"
echo "Sort the output BAM file:                     $sort_bam"
echo "Sort the output BAM file with samtools:       $samtools_sort"
echo "Index the output BAM file:                    $index_bam"
echo "Sample ID (as inferred by the script):        $sampleID"
echo "Output arg for STAR:                          $output_arg"
[[ "$more_args" != "" ]] && echo "Additional args to pass to STAR:              $more_args"
[[ "$annot" != "" ]] && echo "Annotation arg for STAR:                      $annot_arg"
echo "# Listing the input file(s):"
ls -lh "$R1_in"
[[ "$R2_in" != "" ]] && ls -lh "$R2_in"
[[ "$annot" != "" ]] && ls -lh "$annot"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Now creating the output directories...."
mkdir -pv "$outdir"/logs "$outdir"/bam "$starlog_dir"

# Run STAR
echo -e "\n# Now aligning reads with STAR...."
STAR --runThreadN "$threads" \
    --genomeDir "$index_dir" \
    --readFilesIn "$R1_in" "$R2_in" \
    --readFilesCommand zcat \
    --outFilterMultimapNmax $max_map \
    --alignIntronMin $intron_min \
    --alignIntronMax $intron_max \
    --outFileNamePrefix "$outfile_prefix" \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS NM MD \
    --outSAMattrRGline ID:"$sampleID" SM:"$sampleID" \
    $unmapped_arg \
    $output_arg \
    $annot_arg \
    $more_args

#? --twopassMode Basic => Using this following the nf-core RNAseq workflow
#? --outSAMstrandField intronMotif => Using this following the nf-core RNAseq workflow
#? --outSAMattributes NH HI AS NM MD => Using this following the nf-core RNAseq workflow
#TODO
# - Consider using --quantTranscriptomeBan "Singleend", see https://github.com/nf-core/dualrnaseq/blob/master/docs/parameters.md#--quantTranscriptomeBan-Singleend

# Sort with samtools sort
# STAR may fail to sort on some large BAM files, or BAM files with a lot of
# reads mapping to similar positions. In that case, could use samtools sort. 
if [[ "$sort_bam" = true && "$samtools_sort" = "true" ]]; then
    echo -e "\n# Sorting the BAM file with samtools sort..."
    bam_unsorted="$outfile_prefix"Aligned.out.bam
    bam_sorted="$outfile_prefix"Aligned.sortedByCoord.out.bam
    
    samtools sort "$bam_unsorted" > "$bam_sorted"
fi

# Move BAM file
echo -e "\n# Moving the BAM file..."
if [[ "$sort_bam" = true ]]; then
    mv -v "$outfile_prefix"Aligned.sortedByCoord.out.bam "$final_bam"
else
    mv -v "$outfile_prefix"Aligned.out.bam "$final_bam"
fi

# Index BAM file
if [[ "$index_bam" = true ]]; then
    echo -e "\n# Indexing the BAM file..."
    samtools index "$final_bam"
fi

# Organize STAR output
if [ "$output_unmapped" = true ]; then
    echo -e "\n# Moving, renaming and zipping unmapped FASTQ files...."
    for oldpath in "$outfile_prefix"*Unmapped.out.mate*; do
        oldname=$(basename "$oldpath")
        newname=$(echo "$oldname" | sed -E s'/_Unmapped.out.mate([12])/_R\1.fastq.gz/')
        newpath="$unmapped_dir"/"$newname"

        echo "# Moving $oldpath to $newpath..."

        #> The unmapped FASTQ files output by STAR have a weird format with "0:N" for R1 reads (instead of "1:N")
        #> and "1:N" (instead of "2:N") for R2 reads, which Trinity doesn't accept. The code below will fix that:
        [[ "$newpath" = *R1.fastq.gz ]] && sed -E 's/(^@.*) 0:N: (.*)/\1 1:N: \2/' "$oldpath" | gzip -f > "$newpath"
        [[ "$newpath" = *R2.fastq.gz ]] && sed -E 's/(^@.*) 1:N: (.*)/\1 2:N: \2/' "$oldpath" | gzip -f > "$newpath"

        # Remove old file
        rm -v "$oldpath"
    done
fi

# Move STAR log files
echo -e "\n# Moving the STAR log files...."
mv -v "$outfile_prefix"Log*out "$starlog_dir"
echo

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo "========================================================================="
echo -e "# Listing the output BAM file(s):"
ls -lh "$outdir"/bam/"$sampleID"*bam

if [ "$output_unmapped" = true ]; then
    echo "# Listing the FASTQ files with unmapped reads:"
    ls -lh "$unmapped_dir/$sampleID"*fastq.gz
fi

if [ "$count" = true ]; then
    echo -e "\n# Listing the output gene count file:"
    ls -lh "$outdir"/"$sampleID"*ReadsPerGene.out.tab
fi

echo
echo "# STAR version used:"
Print_version | tee "$outdir"/logs/version.txt
echo
[[ "$slurm" = true ]] && Resource_usage
echo
echo "# Done with script"
date
