#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=genemark
#SBATCH --output=slurm-genemark-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               Run GeneMark-ST to predict ORFs in transcripts"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -o <output FASTA> --model <model name> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input FASTA file with transcripts"
    echo "  -o/--outfile    <dir>   Output FASTA file with proteins"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --min_len       <int>   Min protein length (nr of codons)           [default: 100]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to GeneMark-ST"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for CodAn and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/trinity/trinity.fasta -o results/genemark/trinity.faa"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: /fs/ess/PAS0471/jelmer/software/genemark-st/README.GeneMarkS-T"
    echo
}

# Load software
Load_software() {
    GMST=/fs/ess/PAS0471/jelmer/software/genemark-st/gmst.pl
    GMHMMP=/fs/ess/PAS0471/jelmer/software/genemark-st/gmhmmp
    MODEL=/fs/ess/PAS0471/jelmer/software/genemark-st/MetaGeneMark_v1.mod
}
Load_seqkit() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    set -u
}

# Print help for the focal program
Print_help_program() {
    Load_software
    echo "# Help for gm-st:"
    $GMST

    echo -e "\n================================================================="
    echo "# Help for gmHMMp:"
    $GMHMMP
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
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
min_len=100

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
infile=""
outfile=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outfile )    shift && outfile=$1 ;;
        --min_len )         shift && min_len=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
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

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ "$infile" = "" ]] && Die "Please specify an input file with -i/--infile" "$all_args"
[[ "$outfile" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && Die "Input file $infile does not exist"

# Make paths absolute
[[ ! "$infile" =~ ^/ ]] && infile="$PWD"/"$infile"
[[ ! "$outfile" =~ ^/ ]] && outfile="$PWD"/"$outfile"

# Infer the output dir
outdir=$(dirname "$outfile")

# Output file before length-filtering
file_ext=$(basename "$outfile" | sed -E 's/.*(.fasta|.fa|.faa)$/\1/')
outfile_id=$(basename "$outfile" "$file_ext")
outfile_all="$outdir"/"$outfile_id"_all"$file_ext"
infile_id=$(basename "$infile")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT GENEMARK_ST.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $infile"
echo "Output file:                      $outfile"
echo "Min. protein length:              $min_len"
[[ $more_args != "" ]] && echo "Other arguments for Genemark:        $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
${e}mkdir -pv "$outdir"/logs

# Move into the output dir
cd "$outdir" || exit 1

# Run gm-st
echo -e "\n# Running CodAn..."
${e}Time "$GMST" \
    --faa \
    --fnn \
    $more_args \
    "$infile"

#! Threads?

echo -e "\n# Running gmHMMp..."
${e}Time "$GMHMMP" \
    -p 0 \
    -b \
    -m "$MODEL" \
    -A "$infile_id".faa \
    -D "$infile_id".fnn \
    -o "$infile_id".lst \
    "$infile"

echo -e "\n# Renaming the protein FASTA file..."
cp -v "$infile".faa "$outfile_all"

# Remove '*'s at the end - or some programs (e.g. InterProScan) will complain
#echo -e "\n# Removing '*' characters from the proteins..."
#sed -i 's/*$//' "$outfile_all"

# Filter by length
echo "=========================================================================="
echo -e "\n# Filtering the proteins by length with seqkit..."
Load_seqkit
${e}Time seqkit seq \
    --remove-gaps \
    --min-len $min_len \
    "$outfile_all" > "$outfile"

echo "Number of proteins before filtering by length:   $(grep -c "^>" "$outfile_all")"
echo "Number of proteins after filtering by length:    $(grep -c "^>" "$outfile")"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing the final output file:"
    ls -lh "$outfile"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
