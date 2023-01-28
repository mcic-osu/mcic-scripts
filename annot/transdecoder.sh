#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=transdecoder
#SBATCH --output=slurm-transdecoder-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               Run Transdecoder to predict ORFs in transcripts"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -o <output FASTA> --model <model name> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input nucleotide FASTA file with transcripts"
    echo "  -o/--out_aa     <dir>   Output protein FASTA file with proteins"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -O/--out_nuc    <file>  Output nucleotide FASTA - only transcripts with a coding frame"
    echo "  --min_len       <int>   Min protein length (nr of codons)           [default: 100]"
    echo ""
    echo "  --remove_p              Remove the '.p' additions to protein IDs    [default: false]"
    echo "                          Compatible with Evigene output"
    echo "  --keep_stopcod          Keep stop codons '*'s in the sequence       [default: remove]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Transdecoder"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for both the 'gmst' and 'gmhmmp' Genemark utilities, and exit"
    echo "  -v/--version            Print the version of Transdecoder and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/trinity/trinity.fasta -o results/transdecoder/trinity.faa"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/TransDecoder/TransDecoder/wiki"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/transdecoder-5.5.0
    set -u
}
Load_seqkit() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    TransDecoder.LongOrfs --version
    TransDecoder.Predict --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    echo "# Help for TransDecoder.LongOrfs:"
    TransDecoder.LongOrfs
    echo -e "\n================================================================="
    echo "# Help for TransDecoder.Predict:"
    TransDecoder.Predict
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
keep_stopcod=false
remove_p=false

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
infile=""
out_aa=""
out_nuc=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --out_aa )     shift && out_aa=$1 ;;
        -O | --out_nuc )    shift && out_nuc=$1 ;;
        --remove_p )        remove_p=true ;;
        --keep_stopcod )    keep_stopcod=true ;;
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

# Load software
[[ "$dryrun" = false ]] && Load_software

# Check input
[[ "$infile" = "" ]] && Die "Please specify an input file with -i/--infile" "$all_args"
[[ "$out_aa" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && Die "Input file $infile does not exist"

# Make paths absolute
[[ ! "$infile" =~ ^/ ]] && infile=$(realpath "$infile")
[[ ! "$out_aa" =~ ^/ ]] && out_aa="$PWD"/"$out_aa"
[[ ! "$out_nuc" =~ ^/ ]] && out_nuc="$PWD"/"$out_nuc"

# Infer the output dir
outdir=$(dirname "$out_aa")
outdir_nuc=$(dirname "$out_nuc")

# Output file before length-filtering
file_ext=$(basename "$out_aa" | sed -E 's/.*(.fasta|.fa|.faa)$/\1/')
out_aa_id=$(basename "$out_aa" "$file_ext")
out_aa_prelenfilt="$outdir"/"$out_aa_id"_all"$file_ext"
infile_name=$(basename "$infile")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TRANSDECODER.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input nucletotide transcripts:    $infile"
echo "Output protein transcripts:       $out_aa"
[[ "$out_nuc" != "" ]] && echo "Output nucleotide transcripts:    $out_nuc"
echo "Min. protein length:              $min_len"
echo "Keep stop codon symbols:          $keep_stopcod"
echo "Remove '.p' addition to IDs:      $remove_p"
[[ $more_args != "" ]] && echo "Other arguments for Transdecoder: $more_args"
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
${e}mkdir -pv "$outdir"/logs "$outdir_nuc"

# Move into the output dir
cd "$outdir" || exit 1

# Run Transdecoder LongOrfs
echo -e "\n# Running TransDecoder.LongOrfs..."
${e}Time TransDecoder.LongOrfs \
    -m $min_len \
    $more_args \
    -t "$infile"

# Run Transdecoder Predict
echo -e "\n# Running TransDecoder.Predict..."
${e}Time TransDecoder.Predict \
    -t "$infile" \
    --single_best_only

echo -e "\n====================================================================="
# Remove '*'s at the end - or some programs (e.g. InterProScan) will complain
if [[ "$keep_stopcod" = false ]]; then
    echo -e "\n# Removing '*' characters from the proteins..."
    sed 's/*$//' "$infile_name".transdecoder.pep > "$out_aa_prelenfilt"
else
    cp -v "$infile_name".transdecoder.pep "$out_aa_prelenfilt"
fi

# Remove '.p' additions to transcript/protein IDs
if [[ "$remove_p" = true ]]; then
    echo -e "\n# Removing '.p' additions to the protein IDs"
    sed -i -E 's/(t[0-9])\.p[0-9]+ /\1 /' "$out_aa_prelenfilt"
fi

# Filter by length
echo -e "\n# Filtering the proteins by length with seqkit..."
Load_seqkit
${e}Time seqkit seq \
    --remove-gaps \
    --min-len $min_len \
    "$out_aa_prelenfilt" > "$out_aa"

# Filter the nucleotide FASTA to only keep transcripts with a frame
seqkit grep \
    -f <(grep "^>" "$out_aa" | awk '{print $1}' | sed -E 's/>//') \
    <(awk '{print $1}' "$infile") \
    > "$out_nuc"

# Report
echo
echo "Statistics:"
echo "Nr of sequences in the input file:                    $(grep -c "^>" "$infile")"
echo "Nr of proteins with frame, before length-filter:      $(grep -c "^>" "$out_aa_prelenfilt")"
echo "Final nr of proteins (length-filtered, too):          $(grep -c "^>" "$out_aa")"
echo "Final nr of transcripts (should match nr. proteins):  $(grep -c "^>" "$out_nuc")"

# Remove intermediate files
echo -e "\n# Removing intermediate files..."
[[ -s "$out_aa" ]] && rm -v "$out_aa_prelenfilt"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee logs/version.txt
    echo -e "\n# Listing the final output files:"
    ls -lh "$out_aa" "$out_nuc"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
