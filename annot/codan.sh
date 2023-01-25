#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=codan
#SBATCH --output=slurm-codan-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               Run CodAn to predict ORFs in transcripts"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -o <output FASTA> --model <model name> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input FASTA file with transcripts"
    echo "  -o/--outfile    <dir>   Output FASTA file with proteins"
    echo "  --model         <str>   CodAn model"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to CodAn"
    echo "                          One of: FUNGI_full, FUNGI_partial, INV_full, INV_partial, PLANTS_full, PLANTS_partial, VERT_full, VERT_partial"
    echo "                          See https://github.com/pedronachtigall/CodAn/tree/master/models"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for CodAn and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/trinity/trinity.fasta -o results/codan/trinity.faa --model PLANTS_partial"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/pedronachtigall/CodAn"
    echo "  - Docs: https://github.com/pedronachtigall/CodAn/tree/master/tutorial"
    echo "  - Paper: https://academic.oup.com/bib/article/22/3/bbaa045/5847603"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/codan-1.2
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
    Load_software
    codan.py --help
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
        --model )           shift && model=$1 ;;
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

# Infer the output dir
outdir=$(dirname "$outfile")
codan_repo_dir="$outdir"/CodAn
model_dir="$outdir"/"$model"

# Output file before length-filtering
file_ext=$(basename "$outfile" | sed -E 's/.*(.fasta|.fa|.faa)$/\1/')
file_id=$(basename "$outfile" "$file_ext")
outfile_all="$outdir"/"$file_id"_all"$file_ext"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT CODAN.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $infile"
echo "Output file:                      $outfile"
echo "CodAn model:                      $model"
[[ $more_args != "" ]] && echo "Other arguments for CodAn:        $more_args"
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

# Download the models
echo -e "\n# Downloading the CodAn repo and unzipping the desired model..."
[[ ! -d "$codan_repo_dir" ]] && git clone https://github.com/pedronachtigall/CodAn "$codan_repo_dir"
[[ ! -d "$model_dir" ]] && unzip -q -d "$model_dir" "$codan_repo_dir"/models/"$model".zip

# Run CodAn
echo -e "\n# Running CodAn..."
${e}Time codan.py \
    --transcripts="$infile" \
    --model="$model_dir"/"$model" \
    --cpu="$threads" \
    --output="$outdir" \
    $more_args

# Translate sequences #! Note - script assumes a 'partial' model
echo "=========================================================================="
echo -e "\n# Translating the ORF sequences..."
${e}Time TranslatePartial.py "$outdir"/ORF_sequences.fasta "$outfile_all"

# Remove '*'s at the end - or some programs (e.g. InterProScan) will complain
sed -i 's/*$//' "$outfile_all"

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
    ls -lh "$PWD"/"$outfile"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
