#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=star
#SBATCH --output=slurm-star-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Align RNAseq reads to a STAR genome/transcriptome index with STAR
NOTE: STAR is run with several non-default settings, check this script's code for details."
SCRIPT_VERSION="2024-09-25"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=STAR
TOOL_NAME=STAR
TOOL_DOCS="https://github.com/alexdobin/STAR, https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/project/PAS0471/jelmer/conda/star
container_path=/fs/project/PAS0471/containers/depot.galaxyproject.org-singularity-mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2-1df389393721fc66f3fd8778ad938ac711951107-0.img
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true

# Constants - tool parameters
# SEE THE STAR COMMAND BELOW FOR SEVERAL HARDCODED PARAMETERS

# Defaults - tool parameters
quantmode_opt="--quantMode TranscriptomeSAM"    # Will output 'Aligned.toTranscriptome.out.bam' for usage with Salmon 
max_map=20                                      # If this nr is exceeded, read is considered unmapped
sort_bam=samtools
index_bam=false
output_unmapped=false && unmapped_opt=
single_end=false

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage:"
    echo "      sbatch $0 -i data/fastq/S01_R1.fastq.gz -o results/star -r refdata/star_index"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -r/--index_dir      <dir>   Input STAR reference genome index dir (First create index with 'mcic-scripts/rnaseq/star_index.sh')"
    echo "  -o/--outdir         <dir>   BAM output dir"
    echo "  --annot             <file>  Ref. annotation file (GFF/GTF - GTF preferred)"
    echo "                                NOTE: If you don't have or want to use an annotation file,"
    echo "                                omit this this option *and* use the option '--no_transcriptome'."
    echo "  To specify the input reads, use one of the following to options:"
    echo "    -i/--R1           <file>  Input gzipped (R1) FASTQ file path (If R1, name of R2 will be inferred unless using '--single_end')"
    echo "    --fofn            <file>  A File of File Names (FOFN), with one line per input file (not meant for more than 1 sample!)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --no_transcriptome          Don't use '--quantMode TranscriptomeSAM'    [default: use this option, so Salmon can use the output]"
    echo "  --output_unmapped           Output unmapped reads back as FASTQ file    [default: don't output]"
    echo "  --R2                <file>  Input R2 FASTQ file (in case of non-standard naming) [default: infer R2 file name]"
    echo "  --single_end                FASTQ files are single-end                  [default: $single_end]"
    echo "  --sort              <str>   One of 'false', 'star', or 'samtools'       [default: $sort_bam]"
    echo "                              (= Don't sort the BAM file, or sort with STAR or Samtools)"
    echo "  --index_bam                 Index the output BAM file with samtools     [default: $index_bam]"
    echo "  --max_map           <int>   Max. nr. of locations a read can map to     [default: $max_map]"
    echo "  --intron_min        <int>   Min. intron size                            [default: not specified => STAR default]"
    echo "  --intron_max        <int>   Max. intron size                            [default: not specified => STAR default]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you'll have to provide one in order to run the script with a container.)"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --no_strict                 Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
}

# Function to source the script with Bash functions
source_function_script() {
    # Determine the location of this script, and based on that, the function script
    if [[ "$IS_SLURM" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/"$(basename "$FUNCTION_SCRIPT_URL")")
    # Download the function script if needed, then source it
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        function_script=$(basename "$FUNCTION_SCRIPT_URL")
        wget "$FUNCTION_SCRIPT_URL" -O "$function_script"
    fi
    source "$function_script"
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
R1_in= && R2_in= && fofn=
declare -a infiles
index_dir=
annot= && annot_tags= && annot_opt=
intron_min= && intron_min_opt=
intron_max= && intron_max_opt=
outdir=
opts=
version_only=false
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )         shift && outdir=$1 ;;
        -i | --R1 )             shift && R1_in=$1 ;;
        --R2 )                  shift && R2_in=$1 ;;
        --fofn )                shift && fofn=$1 ;;
        -r | --index_dir )      shift && index_dir=$1 ;;
        -a | --annot )          shift && annot=$1 ;;
        --max_map )             shift && max_map=$1 ;;
        --intron_min )          shift && intron_min=$1 ;;
        --intron_max )          shift && intron_max=$1 ;;
        --output_unmapped )     output_unmapped=true ;;
        --no_transcriptome )    quantmode_opt= ;;
        --index_bam )           index_bam=true ;;
        --sort )                shift && sort_bam=$1 ;;
        --single_end )          single_end=true ;;
        --opts )                shift && opts=$1 ;;
        --env )                 shift && env=$1 ;;
        --no_strict )           strict_bash=false ;;
        --dl_container )        dl_container=true ;;
        --container_dir )       shift && container_dir=$1 ;;
        --container_url )       shift && container_url=$1 && dl_container=true ;;
        -h | --help )           script_help; exit 0 ;;
        -v )                    script_version; exit 0 ;;
        --version )             version_only=true ;;
        * )                     die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$R1_in" ]] && die "No input R1 FASTQ file specified, do so with -i/--R1" "$all_opts"
[[ -z "$index_dir" ]] && die "No index dir specified, do so with -r/--index_dir" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$R1_in" ]] && die "Input file $R1_in does not exist"
[[ ! -d "$index_dir" ]] && die "Index dir $index_dir does not exist"
[[ -n "$annot" && ! -f "$annot" ]] && die "Input annotation file (-a) $annot does not exist"
[[ "$sort_bam" != "star" && "$sort_bam" != "samtools" && "$sort_bam" != "false" ]] && die "--sort should be 'false', 'star', or 'samtools', but is $sort_bam"
[[ -n "$quantmode_opt" && -z "$annot" ]] && die "You have not provided an annotation file (which you can do with '--annot'). If you don't want to use an annotation file, you MUST use the '--no_transcriptome' option"

# Input files via FOFN
if  [[ -n "$fofn" ]]; then
    mapfile -t infiles <"$fofn"
    R1_in=${infiles[0]}
    [[ ${#infiles[@]} -eq 2 ]] && R2_in=${infiles[1]}
    [[ ${#infiles[@]} -gt 2 ]] && die "FOFN should contain 1 or 2 filenames, not ${#infiles[@]}"
fi

# Determine R2 file, output prefix, etc
if [[ -n "$R1_in" ]]; then
    R1_basename=$(basename "$R1_in" | sed -E 's/.fa?s?t?q.gz//')
    
    if [[ "$single_end" == false ]]; then
        R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
        sampleID=${R1_basename/"$R1_suffix"/}

        # Determine name of R2 file
        if [[ -z "$R2_in" ]]; then
            R2_suffix=${R1_suffix/1/2}
            R2_in=${R1_in/$R1_suffix/$R2_suffix}
        fi

        [[ ! -f "$R1_in" ]] && die "Input file R1_in $R1_in does not exist"
        [[ ! -f "$R2_in" ]] && die "Input file R2_in $R2_in does not exist"
        [[ "$R1_in" == "$R2_in" ]] && die "Input file R1 is the same as R2"
    else
        sampleID="$R1_basename"
    fi
fi

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
outfile_prefix="$outdir/$sampleID"_
starlog_dir="$outdir"/star_logs
final_bam="$outdir"/bam/"$sampleID".bam
map2trans_bam="$outdir"/bam/"$sampleID"_map2trans.bam
mkdir -p "$outdir"/bam "$starlog_dir"

# If a GFF/GTF file is provided, build the appropriate argument for STAR
if [[ -n "$annot" ]]; then
    annot_opt="--sjdbGTFfile $annot"
    
    if [[ "$annot" =~ .*\.gff3? ]]; then
        annot_tags="--sjdbGTFtagExonParentTranscript Parent"
    elif [[ "$annot" =~ .*\.gtf ]]; then
        annot_tags=
        #These are the defaults, so no need to specify
        #annot_tags="--sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFtagExonParentGene gene_id"
    else
        die "Unknown annotation file format"
    fi
fi

# Sorted output or not
if [[ "$sort_bam" == "star" ]]; then
    output_opt="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 100"
else
    output_opt="--outSAMtype BAM Unsorted"
fi

# Output unmapped reads in a FASTQ file
if [[ "$output_unmapped" == true ]]; then
    unmapped_opt="--outReadsUnmapped Fastx"
    unmapped_dir="$outdir"/unmapped
    mkdir -p "$unmapped_dir"
fi

# Other options
[[ -n "$intron_min" ]] && intron_min_opt="--alignIntronMin $intron_min"
[[ -n "$intron_max" ]] && intron_max_opt="--alignIntronMax $intron_max"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:            $all_opts"
echo "Output BAM dir:                               $outdir"
echo "Input R1 FASTQ file:                          $R1_in"
[[ -n "$R2_in" ]] && echo "Input R2 FASTQ file:                          $R2_in"
echo "Are FASTQ reads single-end?                   $single_end"
echo "Input STAR genome index dir:                  $index_dir"
[[ -n "$annot" ]] && echo "Input annotation file:                        $annot"
echo "Output unmapped reads as FASTQ:               $output_unmapped"
echo "Max nr of alignments for a read:              $max_map"
[[ -n "$intron_min" ]] && echo "Minimum intron size:                          $intron_min"
[[ -n "$intron_max" ]] && echo "Maximum intron size (0 => STAR default):      $intron_max"
echo "Sort the output BAM file:                     $sort_bam"
echo "Index the output BAM file:                    $index_bam"
echo "Sample ID (as inferred by the script):        $sampleID"
echo "Output arg for STAR:                          $output_opt"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
log_time "Listing the input file(s):"
ls -lhd "$index_dir"
ls -lh "$R1_in"
[[ -n "$R2_in" ]] && ls -lh "$R2_in"
[[ -n "$annot" ]] && ls -lh "$annot"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Align with STAR
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    --genomeDir "$index_dir" \
    --readFilesIn "$R1_in" "$R2_in" \
    --outFilterMultimapNmax $max_map \
    --outFileNamePrefix "$outfile_prefix" \
    --runThreadN "$threads" \
    --outSAMattrRGline ID:"$sampleID" SM:"$sampleID" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS NM MD \
    --runRNGseed 0 \
    --alignSJDBoverhangMin 1 \
    --quantTranscriptomeSAMoutput "BanSingleEnd" \
    $intron_min_opt \
    $intron_max_opt \
    $annot_opt \
    $annot_tags \
    $quantmode_opt \
    $unmapped_opt \
    $output_opt \
    $opts

#? --twopassMode Basic => Using this following the nf-core RNAseq workflow
#? --outSAMstrandField intronMotif => Using this following the nf-core RNAseq workflow
#? --outSAMattributes NH HI AS NM MD => Using this following the nf-core RNAseq workflow
#? --quantTranscriptomeSAMoutput "BanSingleEnd" => Using this following the nf-core RNAseq workflow
#? --alignSJDBoverhangMin 1 => Using this following the nf-core RNAseq workflow
#? --runRNGseed 0 => Using this following the nf-core RNAseq workflow

# Sort BAM with samtools sort
#   (STAR may fail to sort some large BAM files, or BAM files with a lot of
#    reads mapping to similar positions. In such cases, use 'samtools sort' instead.)
if [[ "$sort_bam" == "samtools" ]]; then
    log_time "Sorting the main BAM file with samtools sort..."
    bam_unsorted="$outfile_prefix"Aligned.out.bam
    bam_sorted="$outfile_prefix"Aligned.sortedByCoord.out.bam
    runstats samtools sort -o "$bam_sorted" "$bam_unsorted"
    [[ -s "$bam_sorted" ]] && rm -v "$bam_unsorted"
fi

# ==============================================================================
#                           ORGANIZE THE OUTPUT
# ==============================================================================
# Move the output BAM file(s)
log_time "Moving the output BAM file(s)..."
if [[ "$sort_bam" != "false" ]]; then
    mv -v "$outfile_prefix"Aligned.sortedByCoord.out.bam "$final_bam"
else
    mv -v "$outfile_prefix"Aligned.out.bam "$final_bam"
fi
if [[ -n "$quantmode_opt" ]]; then
    mv -v "$outfile_prefix"Aligned.toTranscriptome.out.bam "$map2trans_bam"
fi

# Index the output BAM file(s)
if [[ "$index_bam" == true ]]; then
    log_time "Indexing the output BAM file(s)..."
    runstats samtools index "$final_bam"
    [[ -n "$quantmode_opt" ]] && runstats samtools index "$map2trans_bam"
fi

# Organize unmapped FASTQ files
if [[ "$output_unmapped" == true ]]; then
    log_time "Moving, renaming and zipping unmapped FASTQ files...."
    for oldpath in "$outfile_prefix"*Unmapped.out.mate*; do
        oldname=$(basename "$oldpath")
        newname=$(echo "$oldname" | sed -E s'/_Unmapped.out.mate([12])/_R\1.fastq.gz/')
        newpath="$unmapped_dir"/"$newname"

        log_time "Fixing FASTQ format for $oldpath and outputting $newpath..."
        #> The unmapped FASTQ files output by STAR have a weird format with "0:N" for R1 reads (instead of "1:N")
        #> and "1:N" (instead of "2:N") for R2 reads, which Trinity doesn't accept. The code below will fix that:
        [[ "$newpath" = *R1.fastq.gz ]] && sed -E 's/(^@.*) 0:N: (.*)/\1 1:N: \2/' "$oldpath" | gzip -f > "$newpath"
        [[ "$newpath" = *R2.fastq.gz ]] && sed -E 's/(^@.*) 1:N: (.*)/\1 2:N: \2/' "$oldpath" | gzip -f > "$newpath"
        rm  "$oldpath" # Remove old file
    done
fi

# Move STAR log files
log_time "Moving the STAR log files..."
mv -v "$outfile_prefix"Log*out "$starlog_dir"

# ==============================================================================
#                                   REPORT
# ==============================================================================
# Show alignment % lines from STAR log
log_time "Showing alignment rate lines from $starlog_dir/${sampleID}_Log.final.out..."
grep "Uniquely mapped reads %" "$starlog_dir"/"$sampleID"_Log.final.out
grep "% of reads mapped to multiple loci" "$starlog_dir"/"$sampleID"_Log.final.out

# List output files
log_time "Listing the output BAM file(s):"
ls -lh "$final_bam"
[[ -n "$quantmode_opt" ]] && ls -lh "$map2trans_bam"
if [[ "$output_unmapped" == true ]]; then
    log_time "Listing the FASTQ files with unmapped reads:"
    ls -lh "$unmapped_dir/$sampleID"*fastq.gz
fi

# Final reporting
final_reporting "$LOG_DIR"
