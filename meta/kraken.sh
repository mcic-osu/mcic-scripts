#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=kraken
#SBATCH --output=slurm-kraken-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Kraken2 to assign taxonomy to sequences in a FASTA/FASTQ/pair of FASTQ file(s)"
SCRIPT_VERSION="2023-12-03"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=kraken2
TOOL_NAME=Kraken2
TOOL_DOCS=https://github.com/DerrickWood/kraken2
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/kraken2
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                  # When true, just print tool & script version info and exit

# Defaults - tool(-related) parameters
min_conf=0.15                       # Following https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000949
min_hitgroups=3                     # Following https://www.nature.com/articles/s41596-022-00738-y
write_classif=false                 # Don't write output file(s) with classified reads
write_unclassif=false               # Don't write output file(s) with unclassified reads
use_ram=true                        # Load Kraken db into memory
single_end=false                    # Assume paired-end reads, if the input is FASTQ

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/S1_R1.fastq.gz -o results/kraken --db /fs/project/PAS0471/jelmer/refdata/kraken/std"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile             <file>  Input sequence file (FASTA, single-end FASTQ, or R1 from paired-end FASTQ)"
    echo "                                    - If an R1 paired-end FASTQ file is provided, the name of the R2 file will be inferred"
    echo "                                    - FASTA files should be unzipped; FASTQ files should be gzipped"
    echo "  -o/--outdir             <dir>   Output dir (will be created if needed)"
    echo "  -d/--db                 <dir>   Directory with an existing Kraken database"
    echo "                                    - A few databases are available at: /fs/ess/PAS0471/jelmer/refdata/kraken"
    echo "                                    - Kraken databases can be downloaded from: https://benlangmead.github.io/aws-indexes/k2"
    echo "                                    - Finally, you can use 'mcic-scripts/meta/kraken-build.sh' to create a database"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --confidence            <num>   Confidence required for assignment: number between 0 and 1            [default: $min_conf]"
    echo "  --minimum-hit-groups    <int>   Minimum nr of hit groups required an assignment                       [default: $min_hitgroups]"
    echo "  --memory-mapping                Don't load the full database into RAM memory                          [default: load into memory]"
    echo "                                    - Considerably lower, but useful/needed with very large databases"
    echo "  --classified-out                Write 'classified' sequences to file in '<outdir>/classified' dir     [default: don't write]"
    echo "                                    NOTE: only implemented for PE FASTQ files"
    echo "  --unclassified-out              Write 'unclassified' sequences to file in '<outdir>/unclassified' dir [default: don't write]"
    echo "                                    NOTE: only implemented for PE FASTQ files"
    echo "  --single-end                    FASTQ files are single-end                                            [default: paired-end]"
    echo "  --more_opts             <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env                   <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                    (NOTE: If no default '--container_url' is listed below,"
    echo "                                    you'll have to provide one in order to run the script with a container.)"
    echo "  --conda_env             <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url         <str>   URL to download the container from      [default: $container_url]"
    echo "                                  A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir         <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container                  Force a redownload of the container     [default: $dl_container]"
    echo "  --no_strict                     Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
    echo "  -h/--help                       Print this help message and exit"
    echo "  -v                              Print the version of this script and exit"
    echo "  --version                       Print the version of $TOOL_NAME and exit"
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
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script"
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
infile=
outdir=
db=
conf_opt=
hitgroup_opt=
classif_opt=
unclassif_opt=
mem_opt=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -d | --db | --db-dir )  shift && db=$1 ;;
        --confidence )          shift && min_conf=$1 ;;
        --minimum-hit-groups )  shift && min_hitgroups=$1 ;;
        --single-end )          single_end=true ;;
        --classified-out )      write_classif=true ;;
        --unclassified-out )    write_unclassif=true ;;   
        --memory-mapping )      use_ram=false ;;
        --more_opts )           shift && more_opts=$1 ;;
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
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ -z "$db" ]] && die "No Kraken database dir specified, do so with -d/--db" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -d "$db" ]] && die "Kraken database dir $db does not exist"

# Other prep based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
[[ "$use_ram" = false ]] && mem_opt="--memory-mapping "
[[ -n "$min_conf" ]] && conf_opt="--confidence $min_conf "
[[ -n "$min_hitgroups" ]] && hitgroup_opt="--minimum-hit-groups $min_hitgroups "

# Make sure input file argument is correct based on file type 
if [[ "$infile" =~ \.fa?s?t?q.gz$ ]]; then
    R1_in="$infile"
    R1_basename=$(basename "$R1_in" | sed -E 's/\.fa?s?t?q\.gz//')
    R1_suffix=$(echo "$R1_basename" | sed -E 's/.*(_R?[12]).*/\1/')
    
    if [[ "$single_end" = false ]]; then
        file_type=pe
        R2_suffix=${R1_suffix/1/2}
        R2_in=${R1_in/$R1_suffix/$R2_suffix}
        sample_id=${R1_basename/"$R1_suffix"/}
        infile_opt="--gzip-compressed --paired $R1_in $R2_in"

        [[ ! -f "$R2_in" ]] && die "R2 file $R2_in does not exist"
        [[ "$R1_in" == "$R2_in" ]] && die "R1 file $R1_in is the same as R2 file $R2_in"

        if [[ "$write_classif" = true ]]; then
            classif_opt="--classified-out $outdir/classified/$sample_id#.fastq "
        fi
        if [[ "$write_unclassif" = true ]]; then
            unclassif_opt="--unclassified-out $outdir/unclassified/$sample_id#.fastq "
        fi
    else
        file_type=se
        sample_id=${R1_basename/"$R1_suffix"/}
        infile_opt="--gzip-compressed $R1_in"

        if [[ "$write_classif" == true ]]; then
            classif_opt="--classified-out $outdir/classified/$sample_id.fastq "
        fi
        if [[ "$write_unclassif" == true ]]; then
            unclassif_opt="--unclassified-out $outdir/unclassified/$sample_id.fastq "
        fi
    
    fi
elif [[ "$infile" =~ \.fn?a?s?t?a$ ]]; then
    file_type=fasta
    infile_basename=$(basename "$infile")
    sample_id=${infile_basename%%.*}
    infile_opt="$infile"

    if [[ "$write_classif" == true ]]; then
        classif_opt="--classified-out $outdir/classified/$sample_id.fa "
    fi
    if [[ "$write_unclassif" == true ]]; then
        unclassif_opt="--unclassified-out $outdir/unclassified/$sample_id.fa "
    fi
else
    die "Unknown input file type"
fi

# Define output text files
outfile_main="$outdir"/"$sample_id".main.txt
outfile_report="$outdir"/"$sample_id".report.txt

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output dir:                               $outdir"
echo "Input file type:                          $file_type"
echo "Kraken database dir:                      $db"
echo "Write classified reads?                   $write_classif"
echo "Write unclassified reads?                 $write_unclassif"
[[ -n $min_conf ]] && echo "Minimum confidence:                       $min_conf"
[[ -n $min_hitgroups ]] && echo "Minimum nr of hit groups:                 $min_hitgroups"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
[[ -n $R1_in ]] && ls -lh "$R1_in"
[[ -n $R2_in ]] && ls -lh "$R2_in"
[[ -n $infile ]] && ls -lh "$infile"
ls -lh "$db"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create output dirs
[[ "$write_classif" = true ]] && mkdir -p "$outdir"/classified
[[ "$write_unclassif" = true ]] && mkdir -p "$outdir"/unclassified
mkdir -p "$outdir"/logs

# Run Kraken
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    ${more_opts}--threads "$threads" \
    ${mem_opt}${hitgroup_opt}${conf_opt}--report-minimizer-data \
    ${classif_opt}--db "$db" \
    ${unclassif_opt}--report "$outfile_report" \
    ${infile_opt}> "$outfile_main"

# Rename and zip FASTQ files - only implemented for PE FASTQ
if  [[ "$file_type" = "pe" ]]; then
    if [[ "$write_classif" = true ]]; then
        echo -e "\n# Zipping FASTQ files with classified reads..."
        mv "$outdir"/classified/"$sample_id"_1.fastq "$outdir"/classified/"$sample_id"_R1.fastq
        gzip -f "$outdir"/classified/"$sample_id"_R1.fastq

        mv "$outdir"/classified/"$sample_id"_2.fastq "$outdir"/classified/"$sample_id"_R2.fastq
        gzip -f "$outdir"/classified/"$sample_id"_R2.fastq
    fi
    if [[ "$write_unclassif" = true ]]; then
        echo -e "\n# Zipping FASTQ files with unclassified reads..."
        mv "$outdir"/unclassified/"$sample_id"_1.fastq "$outdir"/unclassified/"$sample_id"_R1.fastq
        gzip -f "$outdir"/unclassified/"$sample_id"_R1.fastq

        mv "$outdir"/unclassified/"$sample_id"_2.fastq "$outdir"/unclassified/"$sample_id"_R2.fastq
        gzip -f "$outdir"/unclassified/"$sample_id"_R2.fastq
    fi
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
