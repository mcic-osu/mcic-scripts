#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=filter_expr
#SBATCH --output=slurm-filter_expr-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "            FILTER A TRANSCRIPTOME BY EXPRESSION LEVEL"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input assembly> -o <output assembly> --kallisto_dir <dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--in_alltrans_nuc    <file>  Input nucleotide FASTA with a transcriptome assembly"
    echo "  -o/--out_alltrans_nuc   <dir>   Output FASTA assembly"
    echo "  --kallisto_dir          <dir>   Base Kallisto dir with Kallisto output"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --out_1trans_nuc        <dir>   Output FASTA assembly with one transcript per gene"
    echo "  --in_1trans_aa          <file>  Input protein FASTA with 1 transcript per gene"
    echo "  --out_1trans_aa         <file>  Output protein FASTA with 1 transcript per gene"
    echo "  --no_genesum                    Don's sum counts by gene (i.e., across transcripts [default: sum by gene]"
    echo "  --min_tpm               <int>   Min. per-sample TPM threshold                      [default: 1]"
    echo "  --mean_tpm              <int>   Mean TPM threshold                                 [default: 0.1]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for filter_low_expr_transcripts and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/trinity/assembly.fasta -o results/trinity/assembly_filt.fasta --kallisto_dir results/kallisto"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#building-expression-matrices"
    echo
}

# Load software
Load_trinity() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/trinity-2.13.2
    set -u
}
Load_seqkit() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    set -u
}
Load_R() {
    module load R/4.2.1-gnu11.2
    FIND_LOW_EXPR=../mcic-scripts/assembly_trans/find_low_expr.R #!TODO FIX
    #[[ ! -f "$FIND_LOW_EXPR" ]] && #TODO - Download the script
}

# Print help for the focal program
Print_help_program() {
    Load_trinity
    filter_low_expr_transcripts.pl
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
min_tpm=1
mean_tpm=0.1
sum_by_gene=true
gene2trans_arg=none

debug=false
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
in_alltrans_nuc=""
out_alltrans_nuc=""
out_1trans_nuc=""
outdir=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --in_alltrans_nuc )    shift && in_alltrans_nuc=$1 ;;
        -o | --out_alltrans_nuc )   shift && out_alltrans_nuc=$1 ;;
        --out_1trans_nuc )          shift && out_1trans_nuc=$1 ;;
        --in_1trans_aa )            shift && in_1trans_aa=$1 ;;
        --out_1trans_aa )           shift && out_1trans_aa=$1 ;;
        --kallisto_dir )            shift && kallisto_dir=$1 ;;
        --min_tpm )                 shift && min_tpm=$1 ;;
        --mean_tpm )                shift && mean_tpm=$1 ;;
        --no_genesum )              sum_by_gene=false ;;
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

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash script settings
set -euo pipefail

# Load software
Load_trinity

# Check input
[[ "$kallisto_dir" = "" ]] && Die "Please specify an input Kallisto dir with --kallisto_dir" "$all_args"
[[ "$in_alltrans_nuc" = "" ]] && Die "Please specify an input assembly with --in_alltrans_nuc" "$all_args"
[[ "$in_alltrans_nuc" = "" ]] && Die "Please specify an output assembly with -o/--in_alltrans_nuc" "$all_args"
[[ ! -f "$in_alltrans_nuc" ]] && Die "Input file $in_alltrans_nuc does not exist"
[[ ! -d "$kallisto_dir" ]] && Die "Input dir $kallisto_dir does not exist"
[[ "$in_1trans_aa" != "" && ! -f "$in_1trans_aa" ]] && Die "Input file $in_1trans_aa does not exist"

# Get outdir, file ID
outdir=$(dirname "$out_alltrans_nuc")
outdir_onetrans=$(dirname "$out_1trans_nuc")
assembly_basename=$(basename "$in_alltrans_nuc")
assembly_id=${assembly_basename/%.*}

# Sum across transcripts?
if [[ "$sum_by_gene" = true ]]; then
    gene2trans="$outdir"/gene2trans.tsv
    gene2trans_arg="--gene_trans_map $gene2trans"
fi

# The abundance_estimates_to_matrix script will create a matrix file with this name:
count_matrix="$outdir"/"$assembly_id".gene.TPM.not_cross_norm

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT FILTER_EXPR.SH"
date
echo "=========================================================================="
echo "All arguments to this script:                     $all_args"
echo
echo "Input assembly (nucl., all transcripts):          $in_alltrans_nuc"
echo "Output assembly (nucl., all transcripts):         $out_alltrans_nuc"
[[ "$out_1trans_nuc" != "" ]] && echo "Output assembly (nucl., 1 transcript per gene):   $out_1trans_nuc"
[[ "$in_1trans_aa" != "" ]] && echo "Input assembly (protein, 1 transcript per gene):  $in_1trans_aa"
[[ "$out_1trans_aa" != "" ]] && echo "Output assembly (protein, 1 transcript per gene): $out_1trans_aa"
echo
echo "Output count matrix:                              $count_matrix"
echo "Sum by gene (across transcripts)?:                $sum_by_gene"
echo "Min. per-sample TPM threshold:                    $min_tpm"
echo "Mean TPM threshold:                               $mean_tpm"
echo
echo "# Listing the input file(s):"
ls -lh "$in_alltrans_nuc" "$kallisto_dir"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs "$outdir_onetrans"

# Create a gene-to-transcript map
if [[ "$sum_by_gene" = true ]]; then
    echo -e "\n# Creating a gene-to-transcript map..."
    Time paste <(grep "^>" "$in_alltrans_nuc" | sed -E 's/>([^ ]+) .*/\1/' | sed -E 's/t[0-9]+//') \
        <(grep "^>" "$in_alltrans_nuc" | sed -E 's/>([^ ]+) .*/\1/') |
        sort -k1,1 > "$gene2trans"
fi

# Run abundance_estimates_to_matrix.pl
if [[ ! -f "$count_matrix" ]]; then
    echo -e "\n# Running abundance_estimates_to_matrix..."
    Time abundance_estimates_to_matrix.pl \
        --est_method kallisto \
        $gene2trans_arg \
        --out_prefix "$outdir"/"$assembly_id" \
        --name_sample_by_basedir \
        $(find "$kallisto_dir" -name "abundance.tsv")
else
    echo -e "\n# NOTE: Skipping matrix tabulation, file exists..."
    ls -lh "$count_matrix"
fi

# Filter by expression level
echo -e "\n# Running script find_low_expr.R to get IDs of low-expression genes..."
Load_R
Rscript "$FIND_LOW_EXPR" "$count_matrix" "$outdir"/lowTPM_geneIDs.txt $min_tpm $mean_tpm

join "$gene2trans" <(sort "$outdir"/low_TPM_ids.txt) | \
    cut -d " " -f 2 > "$outdir"/lowTPM_transIDs.txt # Get transcript IDs (from gene IDs)

echo -e "\n# Removing low-expression genes from the assembly..."
Load_seqkit
Time seqkit grep -v \
    -f "$outdir"/lowTPM_transIDs.txt \
    <(awk '{print $1}' "$in_alltrans_nuc") \
    > "$out_alltrans_nuc"

# Create an assembly with one transcript per gene
if [[ "$out_1trans_nuc" != "" ]]; then
    echo -e "\n# Creating a transcriptome with 1 transcript per gene..."
    Time seqkit grep --by-name -r -p "t1$" "$out_alltrans_nuc" > "$out_1trans_nuc"
fi

# Also filter the AA FASTA
#? NOTE: Fasta IDs must contain a space + 2nd word, or final Diamond output will have 'Query_1' etc as IDs!)
echo -e "\n# Filtering the protein FASTA file to remove low-expression genes ..."
seqkit grep \
    -f <(grep "^>" "$out_1trans_nuc" | awk '{print $1}' | sed -E 's/>//') \
    <(awk '{print $1}' "$in_1trans_aa") |
    sed '/^>/s/$/ frame/' \
    >"$out_1trans_aa" 

# Count the number of transcripts
echo
echo "# In- and output statistics:"
[[ "$in_1trans_aa" != "" ]] && echo "Count for input - aa - 1trans:                $(grep -c "^>" "$in_1trans_aa")"
echo "Count for input - nuc - all-transcripts:      $(grep -c "^>" "$in_alltrans_nuc")"
echo "Count for output - nuc - all-transcripts:     $(grep -c "^>" "$out_alltrans_nuc")"
[[ "$out_1trans_nuc" != "" ]] && echo "Count for output - nuc - 1trans:              $(grep -c "^>" "$out_1trans_nuc")"
[[ "$out_1trans_aa" != "" ]] && echo "Count for output - aa - 1trans:               $(grep -c "^>" "$out_1trans_aa")"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "# Listing the output assembly:"
ls -lh "$(realpath "$out_alltrans_nuc")"
[[ "$out_1trans_nuc" != "" ]] && ls -lh "$(realpath "$out_1trans_nuc")"
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
echo
