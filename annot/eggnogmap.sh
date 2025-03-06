#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=eggnogmap
#SBATCH --output=slurm-eggnogmap-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run EggNOGmapper to functionally annotate proteins"
readonly MODULE=miniconda3/23.3.1-py310
readonly CONDA=/fs/ess/PAS0471/jelmer/conda/eggnogmapper
readonly SCRIPT_VERSION="2023-06-30"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY=emapper.py
readonly TOOL_NAME=EggNOGmapper
readonly TOOL_DOCS=https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.11
readonly VERSION_COMMAND="$TOOL_BINARY --version | grep emapper"

# Constants - parameters
INPUT_TYPE=proteins             # Assume a protein (rather than nucleotide) FASTA

# Parameter defaults
search_method=diamond           # Also the EggNOGmapper default
go_evidence="non-electronic"    # Also the EggNOGmapper default
sensmode="more-sensitive"       # EggNOGmapper default is 'sensitive'
tax_scope="auto"                # Also the EggNOGmapper default

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 -i results/braker/proteins.faa -d data/eggnog -o results/eggnogmapper"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--proteins   <file>  Input protein FASTA file"
    echo "  -d/--db_dir     <dir>   Pre-downloaded eggNOGmapper database dir"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --out_prefix    <str>   Output prefix, e.g. a genome ID             [default: basename of protein file]"
    echo "  --tax_scope     <str>   Taxonomic scope (see eggNOGmapper docs)     [default: 'auto']"
    echo "  --search_method <str>   'diamond', 'mmseqs', or 'hmmer'             [default: 'diamond']"
    echo "  --sensmode      <str>   Diamond sensitivity                         [default: 'more-sensitive']"
    echo "                            One of 'default', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive' or 'ultra-sensitive'"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to EggNOGmapper"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - This script is set up to work with a protein FASTA file;"
    echo "    whereas it is also possible to run EggNOGmapper with a nucleotide FASTA file."
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    local is_slurm=$1

    # Determine the location of this script, and based on that, the function script
    if [[ "$is_slurm" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/bash_functions.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
    source "$function_script"
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
db_dir=
out_prefix=
outdir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && readonly infile=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        -d | --db_dir )     shift && readonly db_dir=$1 ;;
        --out_prefix )      shift && readonly out_prefix=$1 ;;
        --search_method )   shift && readonly search_method=$1 ;;
        --sensmode )        shift && readonly sensmode=$1 ;;
        --tax_scope )       shift && readonly tax_scope=$1 ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version )         load_env "$MODULE" "$CONDA"
                            tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$db_dir" ]] && Die "No EggNOGmapper DB dir specified, do so with -d/--db_dir" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -d "$db_dir" ]] && die "Input EggNOGmapper DB dir $db_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Other variables
if [[ -z "$out_prefix" ]]; then
    prot_basename=$(basename "$infile")
    out_prefix=${prot_basename%.*}
fi

if [[ "$IS_SLURM" == true ]]; then
    scratch_dir="$PFSDIR"
    temp_dir="$PFSDIR"/tmp
    tempdir_arg="--scratch_dir $scratch_dir --temp_dir $temp_dir"
else
    temp_dir="$outdir"/tmp
    tempdir_arg="--temp_dir $temp_dir"
fi

# Get nr of genes in the input
ngenes_in=$(grep -c "^>" "$infile")

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input file:                                   $infile"
echo "Nr of entries (genes) in the input file:      $ngenes_in"
echo "eggNOGmapper database dir:                    $db_dir"
echo "Output dir:                                   $outdir"
echo "Output prefix:                                $out_prefix"
echo "Search method:                                $search_method"
echo "DIAMOND sensitivity:                          $sensmode"
echo "Taxonomic scope:                              $tax_scope"
echo "Scratch and temp dir argument:                $tempdir_arg"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
mkdir -p "$temp_dir"

# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
        -i "$infile" \
        --itype "$INPUT_TYPE" \
        --data_dir "$db_dir" \
        --output_dir "$outdir" \
        --output "$out_prefix" \
        -m "$search_method" \
        --sensmode "$sensmode" \
        --tax_scope "$tax_scope" \
        --go_evidence "$go_evidence" \
        --cpu "$threads" \
        --override \
        $tempdir_arg \
        $more_args

# Report
ngenes_out=$(grep -cv "^#" "$outdir"/*emapper.annotations)
ngenes_descrip=$(grep -v "^#" "$outdir"/*emapper.annotations | cut -f 8 | grep -cv "^-")
echo
echo "Nr of genes/entries in the input:         $ngenes_in"
echo "Nr of genes/entries in the output:        $ngenes_out"
echo "Nr of genes/entries with a description:   $ngenes_descrip"

#? Non-default options for Eggnogmapper
# --go_evidence all => default is to use only non-electronic terms (`non-electronic`), see https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8
# --go_evidence {experimental,non-electronic,all}
#                        Defines what type of GO terms should be used for annotation. experimental = Use only terms inferred from experimental evidence. non-electronic = Use only non-electronically curated terms
# --override => Overwrite existing output files

#? Other options for Eggnogmapper
# --pfam_realign denovo Needs some HMMer server setup
#--list_taxa            List taxa available for --tax_scope/--tax_scope_mode, and exit
#--tax_scope            ....
#--resume               Resumes a previous emapper run, skipping results in existing output files.

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
