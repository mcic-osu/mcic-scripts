#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=eggnogmap
#SBATCH --output=slurm-eggnogmap-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
search_method=diamond           # Also the EggNOGmapper default
go_evidence="non-electronic"    # Also the EggNOGmapper default
sensmode="more-sensitive"       # EggNOGmapper default is 'sensitive'
tax_scope="auto"                # Also the EggNOGmapper default
slurm=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "        Run eggNOGmapper to functionally annotate proteins"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -d <database dir> -o <output dir> [...]"
    echo "  bash $0 -h"
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
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for EggNOGmapper and exit"
    echo "  -v/--version            Print the version of EggNOGmapper and exit"
    echo
    echo "HARDCODED OPTIONS:"
    echo "  - This script is set up to work with a protein FASTA file;"
    echo "    whereas it is also possible to run EggNOGmapper with a nucleotide FASTA file."
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/braker/proteins.faa -d data/eggnog -o results/eggnogmapper"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /users/PAS0471/jelmer/miniconda3/envs/eggnogg-env
    export EGGNOG_DATA_DIR="$db_dir"
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    emapper.py --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    emapper.py --help
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
    echo "Memory (MB per node): $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):      $SLURM_CPUS_PER_TASK"
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
    date
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
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
proteins=""
db_dir=""
out_prefix=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --proteins )   shift && proteins=$1 ;;
        -d | --db_dir )     shift && db_dir=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --out_prefix )      shift && out_prefix=$1 ;;
        --search_method )   shift && search_method=$1 ;;
        --sensmode )        shift && sensmode=$1 ;;
        --tax_scope )       shift && tax_scope=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        * )                 Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
Load_software
Set_threads

# Check input
[[ "$proteins" = "" ]] && Die "Please specify an input file with -i/--proteins" "$all_args"
[[ "$db_dir" = "" ]] && Die "Please specify an input file with -d/--db_dir" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$proteins" ]] && Die "Input file $proteins does not exist"
[[ ! -d "$db_dir" ]] && Die "Input dir $db_dir does not exist"

# Other variables
if [[ "$out_prefix" = "" ]]; then
    prot_basename=$(basename "$proteins")
    out_prefix=${prot_basename%.*}
fi

if [[ "$slurm" = true ]]; then
    scratch_dir="$PFSDIR"
    temp_dir="$PFSDIR"/tmp
    tempdir_arg="--scratch_dir $scratch_dir --temp_dir $temp_dir"
else
    temp_dir="$outdir"/tmp
    tempdir_arg="--temp_dir $temp_dir"
fi

# Get nr of genes in the input
ngenes_in=$(grep -c "^>" "$proteins")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT EGGNOGMAP.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $proteins"
echo "eggNOGmapper database dir:        $db_dir"
echo "Output dir:                       $outdir"
echo "Output prefix:                    $out_prefix"
echo "Search method:                    $search_method"
echo "DIAMOND sensitivity:              $sensmode"
echo "Taxonomic scope:                  $tax_scope"
echo "Scratch and temp dir argument:    $tempdir_arg"
[[ $more_args != "" ]] && echo "Other arguments for EggNOGmapper: $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Nr of genes/entries in the input: $ngenes_in" 
echo "Listing the input file(s):"
ls -lh "$proteins"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs "$temp_dir"

# Run
echo -e "\n# Running eggNOGmapper..."
Time \
    emapper.py \
        -i "$proteins" \
        --itype proteins \
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

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing files in the output dir:"
ls -lhd "$PWD"/"$outdir"/*
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
echo
