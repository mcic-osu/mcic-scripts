#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=filter_alltrans
#SBATCH --output=slurm-filter_alltrans-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              FILTER AN ALL-TRANSCRIPTS TRANSCRIPTOME"
    echo "    TO ONLY KEEP GENES IN A 1-TRANSCRIPT-PER-GENE TRANSCRIPTOME"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -I <input FASTA> -o <output FASTA> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--asm_alltrans_in   <file>   Input FASTA assembly with all transcripts"
    echo "  -I/--asm_1trans        <file>   Input FASTA assembly with a single transcript per gene"
    echo "  -o/--asm_alltrans_out  <file>   Output assembly FASTA (will have all transcripts)"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                         Run the script in debug mode (print all code)"
    echo "  -h / --help                     Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --entap_dir results/entap --asm_1trans results/assembly.fa --asm_alltrans results/assembly_alltrans.fa \ "
    echo "     --asm_out results/entap_proc/assembly.fa --annot_out results/entap_proc/annot.tsv"
    echo
}

# Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqkit
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
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
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
onetrans=""               # Input assembly with 1 transcript per genes
alltrans_in=""            # Input assembly with all transcripts
alltrans_out=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --asm_alltrans_in )    shift && alltrans_in=$1 ;;
        -I | --asm_1trans )         shift && onetrans=$1 ;;
        -o | --asm_alltrans_out )   shift && alltrans_out=$1 ;;
        -v | --version )            Print_version; exit 0 ;;
        -h )                        Print_help; exit 0 ;;
        --help )                    Print_help_program; exit 0;;
        --debug )                   debug=true ;;
        * )                         Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Load software
Load_software

# Bash script settings
set -euo pipefail

# Output files
outdir=$(dirname "$alltrans_out")
gene2trans="$outdir"/gene2trans.tsv                 # Gene-to-transcript map
trans_final="$outdir"/transIDs_final.txt            # List of transcripts

# Check input
[[ "$onetrans" = "" ]] && Die "Please specify an input assembly with --asm_1trans" "$all_args"
[[ "$alltrans_in" = "" ]] && Die "Please specify an input assembly with --asm_alltrans_in" "$all_args"
[[ "$alltrans_out" = "" ]] && Die "Please specify an output assembly with --asm_alltrans_out" "$all_args"

[[ ! -f "$onetrans" ]] && Die "Input file $onetrans does not exist"
[[ ! -f "$alltrans_in" ]] && Die "Input file $alltrans_in does not exist"

# Report
echo
echo "=========================================================================="
echo "                 STARTING SCRIPT FILTER_ALLTRANS.SH"
date
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input assembly - 1 transcript per gene:   $onetrans"
echo "Input assembly - all transcripts:         $alltrans_in"
echo "Output assembly - all transcripts:        $alltrans_out"
echo
echo "Listing the input files:"
ls -lh "$onetrans" "$alltrans_in"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs

echo -e "\n# Subsetting the assembly with all isoforms..."
# Create a gene-to-transcript map
paste <(grep "^>" "$alltrans_in" | sed -E 's/>([^ ]+) .*/\1/' | sed -E 's/t[0-9]+//') \
    <(grep "^>" "$alltrans_in" | sed -E 's/>([^ ]+) .*/\1/') |
    sort -k1,1 > "$gene2trans"

# Create a list of final transcripts to keep
join -t $'\t' "$gene2trans" \
    <(grep "^>" "$onetrans" | awk '{print $1}' | sed -e 's/t1//' -e 's/>//' | sort) |
    cut -f 2 \
    > "$trans_final"

# Subset the FASTA file
seqkit grep -f "$trans_final" "$alltrans_in" > "$alltrans_out"

# Get stats to report
ngene_in=$(grep -c "^>" "$onetrans")
ntrans_in=$(grep -c "^>" "$alltrans_in")
ntrans_out=$(grep -c "^>" "$alltrans_out")

# Report
echo 
echo "Number of input 'genes':          $ngene_in"
echo "Number of input transcripts:      $ntrans_in"
echo "Number of output transcripts:     $ntrans_out"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "\n# Listing the output files:"
ls -lh "$alltrans_out"
echo -e "\n# Done with script"
date
