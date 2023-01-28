#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
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
    echo "  -i/--assembly_in    <file>  Input nucleotide FASTA with a transcriptome assembly"
    echo "  -o/--assembly_out   <dir>   Output FASTA assembly"
    echo "  --kallisto_dir      <dir>   Base Kallisto dir with Kallisto output"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -O/--onetrans_out   <dir>   Output FASTA assembly with one transcript per gene"
    echo "  --no_genesum                Don's sum counts by gene (i.e., across transcripts [default: sum by gene]"
    echo "  --min_tpm           <int>   Min. per-sample TPM threshold                      [default: 1]"
    echo "  --mean_tpm          <int>   Mean TPM threshold                                 [default: 0.1]"
    echo "  --more_args         <str>   Quoted string with additional argument(s) to pass to filter_low_expr_transcripts.pl"
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
assembly_in=""
assembly_out=""
onetrans_out=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly_in )    shift && assembly_in=$1 ;;
        -o | --assembly_out )   shift && assembly_out=$1 ;;
        -O | --onetrans_out )   shift && onetrans_out=$1 ;;
        --kallisto_dir )        shift && kallisto_dir=$1 ;;
        --min_tpm )             shift && min_tpm=$1 ;;
        --mean_tpm )            shift && mean_tpm=$1 ;;
        --no_genesum )          sum_by_gene=false ;;
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

# Bash script settings
set -euo pipefail

# Load software
Load_trinity

# Check input
[[ "$assembly_in" = "" ]] && Die "Please specify an input assembly with -i/--assembly_in" "$all_args"
[[ "$kallisto_dir" = "" ]] && Die "Please specify an input Kallisto dir with --kallisto_dir" "$all_args"
[[ "$assembly_in" = "" ]] && Die "Please specify an output assembly with -o/--assembly_in" "$all_args"
[[ ! -f "$assembly_in" ]] && Die "Input file $assembly_in does not exist"
[[ ! -d "$kallisto_dir" ]] && Die "Input dir $kallisto_dir does not exist"

# Get outdir, file ID
outdir=$(dirname "$assembly_out")
outdir_onetrans=$(dirname "$onetrans_out")
assembly_basename=$(basename "$assembly_in")
assembly_id=${assembly_basename/%.*}
assembly_intermed="$outdir"/"$assembly_id"_intermed.fa

# Sum across transcripts?
if [[ "$sum_by_gene" = true ]]; then
    gene2trans="$outdir"/gene2trans.tsv
    gene2trans_arg="--gene_trans_map $gene2trans"
fi

# The abundance_estimates_to_matrix script will create a matrix file with this name:
count_matrix="$outdir"/"$assembly_id".isoform.TPM.not_cross_norm

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT FILTER_EXPR.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input assembly:                   $assembly_in"
echo "Output assembly:                  $assembly_out"
[[ "$onetrans_out" = "" ]] && echo "Output assembly with 1 transcript per gene: $onetrans_out"
echo "Sum by gene (across transcripts?) $sum_by_gene"
echo "Min. per-sample TPM threshold:    $min_tpm"
echo "Mean TPM threshold:               $mean_tpm"
[[ $more_args != "" ]] && echo "Other arguments for filter_low_expr_transcripts:  $more_args"
echo
echo "Listing the input file(s):"
ls -lh "$assembly_in" "$kallisto_dir"
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
    paste <(grep "^>" "$assembly_in" | sed -E 's/>([^ ]+) .*/\1/' | sed -E 's/t[0-9]+//') \
        <(grep "^>" "$assembly_in" | sed -E 's/>([^ ]+) .*/\1/') |
        sort -k1,1 > "$gene2trans"
fi

# Run abundance_estimates_to_matrix.pl
echo -e "\n# Running abundance_estimates_to_matrix..."
Time abundance_estimates_to_matrix.pl \
    --est_method kallisto \
    $gene2trans_arg \
    --out_prefix "$outdir"/"$assembly_id" \
    --name_sample_by_basedir \
    $(find "$kallisto_dir" -name "abundance.tsv")

# Filter by min. per-sample TPM (expression level)
echo -e "\n# Running filter_low_expr_transcripts to filter by min. per-sample TPM..."
Time filter_low_expr_transcripts.pl \
    --matrix "$count_matrix" \
    --transcripts "$assembly_in" \
    --min_expr_any $min_tpm \
    $more_args \
    > "$assembly_intermed"

# Get list of genes with too-low mean TPM
echo -e "\n# Filtering by overall mean TPM..."
tail -n +2 "$count_matrix" |
    awk '{s=0; for(i=2; i<=NF; i++) s=s+$i; print $1 "\t" s/(NF-1)}' |
    awk -v mean_tpm="$mean_tpm" '$2 < mean_tpm' | cut -f 1 \
    > "$outdir"/low_mean_TPM_ids.txt

# Filter by mean TPM
Load_seqkit
seqkit grep -v -f "$outdir"/low_mean_TPM_ids.txt "$assembly_intermed" > "$assembly_out"

# Get a file with one transcript per genes
if [[ "$onetrans_out" != "" ]]; then
    echo -e "\n# Creating a transcriptome with 1 transcript per gene..."
    seqkit grep "t1 " "$assembly_out" > "$onetrans_out"
fi

# Count the number of transcripts
echo
echo "Nr transcripts in the input assembly:           $(grep -c "^>" "$assembly_in")"
echo "Nr transcripts after min. TPM filtering:        $(grep -c "^>" "$assembly_intermed")"
echo "Nr transcripts after min+mean. TPM filtering:   $(grep -c "^>" "$assembly_out")"
[[ "$onetrans_out" != "" ]] && echo "Nr final genes:                                   $(grep -c "^>" "$onetrans_out")"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "# Listing the output assembly:"
ls -lh "$(realpath "$assembly_out")"
[[ "$onetrans_out" != "" ]] && ls -lh "$(realpath "$onetrans_out")"
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
echo
