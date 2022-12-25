#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --mem=12G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fqsub
#SBATCH --output=slurm-fqsub-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
        echo
    echo "======================================================================"
    echo "                            $0"
    echo "               Script to subsample FASTQ files using seqtk"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <(R1) FASTQ input> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1_in      <file>  Input FASTQ file"
    echo "                             For paired-end reads, provide R1 file and the script will detect the R1 file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --single_end            Sequences are single-end             [default: paired-end]"
    echo "  --n_reads       <int>   Number of reads to select            [default: 100,000]"
    echo "  --prop_reads    <num>   Proportion of reads to select, e.g. '0.1'"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for seqtk and exit"
    echo "  -v/--version            Print the version of seqtk and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/sampleA_R1.fastq.gz -o results/fq_subset"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqtk
}

## Print version
Print_version() {
    set +eo pipefail
    Load_software
    seqtk 2>&1 | head -n 3 | tail -n 1
    set -eo pipefail
}

## Print help for the focal program
Print_help_program() {
    Load_software
    seqtk
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | grep -Ev "ba|ex"
    echo
}

## Print SLURM job requested resources
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

## Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
rand=$RANDOM

## Option defaults
n_reads=100000
single_end=false

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
R1_in=""
outdir=""
prop_reads=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1_in )      shift && R1_in=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --single_end )      single_end=true ;;
        --n_reads )         shift && n_reads=$1 ;;
        --prop_reads )      shift && prop_reads=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
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
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software

## Bash script settings
set -euo pipefail

## Infer the input dir
indir=$(dirname "$R1_in")
R1_out="$outdir"/$(basename "$R1_in")

## FASTQ filename parsing
if [[ "$single_end" = false ]]; then
    file_ext=$(basename "$R1_in" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)/\1/')
    R1_suffix=$(basename "$R1_in" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    R2_out="$outdir"/$(basename "$R2_in")
fi

## Check input - error out if neither n_reads or prop_reads is provided
[[ $n_reads = "" ]] && [[ $prop_reads = "" ]] && \
    echo "Die neither a nr. (-n) or a proportion (-p) of reads is provided" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$R1_in" ]] && Die "Input file $R1_in does not exist"
[[ "$R2_in" != "" ]] && [[ ! -f "$R2_in" ]] && Die "Input file $R2_in does not exist"
[[ "$indir" = "$outdir" ]] && Die "Input and output dirs can't be the same!"
[[ "$R1_in" = "$R2_in" ]] && Die "Name parsing error: input R1 and R2 are the same file!"

## Number of reads in input FASTQ file
n_reads_total=$(zcat "$R1_in" | awk '{ s++ } END{ print s/4 }')

## If prop_reads is given, calculate n_reads
[[ $prop_reads != "" ]] && n_reads=$(python -c "print(int($n_reads_total * $prop_reads))")

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT FQ_SUB.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input (R1) FASTQ file:            $R1_in"
echo "Output dir:                       $outdir"
[[ $prop_reads != "" ]] && echo "Proportion of reads to keep:      $prop_reads"
echo "Number of reads to keep:          $n_reads"
echo "Reads are single-end:             $single_end"
echo
echo "Random seed:                      $rand"
echo "(R1) output file:                 $R1_out"
[[ "$single_end" = false ]] && echo "R2 output file:                   $R2_out"
echo "Total (input) nr of reads:        $n_reads_total"
echo
echo "Listing the input file(s):"
ls -lh "$R1_in"
[[ "$single_end" = false ]] && ls -lh "$R2_in"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="


# =============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

if [[ "$dryrun" = false ]]; then
    
    ## Run seqtk
    Time seqtk sample -s$rand "$R1_in" "$n_reads" | gzip > "$R1_out"
    [[ "$single_end" = false ]] && Time seqtk sample -s$rand "$R2_in" "$n_reads" | gzip > "$R2_out"

    ## Count reads in output
    n_reads_R1_out=$(zcat "$R1_out" | awk '{ s++ } END{ print s/4 }')
    [[ "$single_end" = false ]] && n_reads_R2_out=$(zcat "$R2_out" | awk '{ s++ } END{ print s/4 }')
    echo -e "\n# Number of reads in output file(s):"
    echo "R1 out: $n_reads_R1_out"
    [[ "$single_end" = false ]] && echo "R2 out: $n_reads_R2_out"

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lh "$PWD"/"$R1_out"
    [[ "$single_end" = false ]] && ls -lh "$PWD"/"$R2_out"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
