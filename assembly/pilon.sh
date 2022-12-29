#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --mem=172G
#SBATCH --cpus-per-task=25
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=pilon
#SBATCH --output=slurm-pilon-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  RUN PILON TO POLISH A GENOME"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --genome <FASTA> --bam_dir <dir with BAMs> --out_prefix <output prefix> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --genome        <file>  Genome assembly FASTA file"
    echo "  --bam_dir       <dir>   Directory with BAM files of Illumina reads mapped to the assembly"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  --out_prefix    <str>   Prefix for output files"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --fix           <str>   What to fix: 'snps'/'indels'/'gaps'/'local'/'all'/'bases'   [default: 'bases']"
    echo "                          See the Pilon documentation for details"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Pilon"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Pilon and exit"
    echo "  -v/--version            Print the version of Pilon and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --genome results/cany/assemblyA.fasta --bam_dir results/star --out_prefix assemblyA -o results/pilon"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - https://github.com/broadinstitute/pilon/wiki"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/pilon-1.24
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    pilon --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    pilon --help
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
fix=bases

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
genome=""
bam_dir="" && bam_arg=""
outdir=""
more_args=""


# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --genome )         shift && genome=$1 ;;
        --bam_dir )        shift && bam_dir=$1 ;;
        -o | --outdir )    shift && outdir=$1 ;;
        --out_prefix )     shift && out_prefix=$1 ;;
        --fix )            shift && fix=$1 ;;
        --more_args )      shift && more_args=$1 ;;
        --debug )          debug=true ;;
        --dryrun )         dryrun=true ;;
        -v | --version )   Print_version; exit ;;
        -h )               Print_help; exit 0 ;;
        --help )           Print_help_program; exit 0;;
        * )                Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Bash script settings
set -euo pipefail

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Load software
[[ "$dryrun" = false ]] && Load_software

# Build BAM arg
for bam in "$bam_dir "/*bam; do
    bam_arg="$bam_arg --frags $bam"
done

# Check input
[[ $genome = "" ]] && Die "Please specify an input genome FASTA file with -i/--genome"
[[ $bam_dir  = "" ]] && Die "Please specify an input BAM dir with -I/--bam"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o/--outdir"
[[ $out_prefix = "" ]] && Die "Please specify an output prefix with -p/--out_prefix"
[[ ! -f $genome ]] && Die "Genome FASTA $genome does not exist"
[[ ! -d $bam_dir  ]] && Die "BAM dir $bam_dir  does not exist"

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT PILON.SH"
date
echo "=========================================================================="
echo "All arguments to this script: $all_args"
echo "Input genome assembly FASTA:  $genome"
echo "Input BAM:                    $bam"
echo "Output dir:                   $outdir"
echo "Output prefix:                $out_prefix"
echo "What to fix (--fix):          $fix"
[[ $more_args != "" ]] && echo "Other arguments for Pilon:    $more_args"
echo
echo "# Listing the input genome FASTA:"
ls -lh "$genome"
echo
echo "# BAM file argument:"
echo "$bam_arg"
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
${e}mkdir -pv "$outdir"/logs

echo -e "\n# Now running Pilon..."
${e}Time pilon \
    --genome "$genome" \
    $bam_arg \
    --outdir "$outdir" \
    --output "$out_prefix" \
    --fix "$fix" \
    $more_args

# No --threads: running with v 1.24, got: "--threads argument no longer supported; ignoring!"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo
echo "# Done with script"
date
