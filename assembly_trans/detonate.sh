#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=detonate
#SBATCH --output=slurm-detonate-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "           EVALUATE A TRANSCRIPTOME ASSEMBLY WITH DETONATE"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input assembly> -I <input-FASTQ-dir> -o <output-dir> ...[...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly   <file>  Input transcriptome assembly FASTA file"
    echo "  -I/--fq-dir     <dir>   Input dir with FASTQ files"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --read-length   <int>   FASTQ read length                   [default: 150]"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Detonate"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Detonate and exit"
    echo "  -v/--version            Print the version of Detonata and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly.fa -I data/fastq/ -o results/detonate"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - The script assumes that the FASTQ reads are paired-end"
    echo
    echo "OUTPUT:"
    echo "  - "
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: http://deweylab.biostat.wisc.edu/detonate/vignette.html"
    echo "  - Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4298084/"
    echo
}

## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /users/PAS0471/jelmer/miniconda3/envs/detonate-env
    set -u
}

## Print version
Print_version() {
    Load_software
    rsem-eval-calculate-score --version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    rsem-eval-calculate-score --help
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
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

## Set the number of threads/CPUs
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
## Option defaults
readlen=150

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
assembly=""
fqdir=""
outdir=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly )   shift && assembly=$1 ;;
        -I | --fq-dir )     shift && fqdir=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --read-length )     shift && readlen=$1 ;;
        --more-args )       shift && more_args=$1 ;;
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
Set_threads

## Bash script settings
set -euo pipefail

## Other parameters
R1_LIST=$(echo "$fqdir"/*R1*fastq.gz | sed 's/ /,/g')
R2_LIST=$(echo "$fqdir"/*R2*fastq.gz | sed 's/ /,/g')
file_ext=$(basename "$assembly" | sed -E 's/.*(.fasta|.fa|.fna)$/\1/')
PREFIX="$outdir"/$(basename "$assembly" "$file_ext")

## Check input
[[ "$assembly" = "" ]] && Die "Please specify an input assembly with -i/--assembly"
[[ "$fqdir" = "" ]] && Die "Please specify a FASTQ input dir with -I/--fqdir"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir"
[[ ! -f "$assembly" ]] && Die "Input assembly file $assembly does not exist"
[[ ! -d "$fqdir" ]] && Die "Input FASTQ dir $fqdir does not exist"

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT DETONATE.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input assembly:                   $assembly"
echo "Input FASTQ dir:                  $fqdir"
echo "Output dir:                       $outdir"
echo "Read length:                      $readlen"
[[ $more_args != "" ]] && echo "Other arguments for Detonate:     $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "List of R1 files:                 $R1_LIST"
echo "List of R2 files:                 $R2_LIST"
echo "Listing the input FASTQ file(s):"
ls -lh "$fqdir"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

echo -e "\nRunning rsem-eval-calculate-score..."
${e}Time rsem-eval-calculate-score \
    --paired-end \
    "$R1_LIST" "$R2_LIST" \
    "$assembly" \
    "$PREFIX" \
    "$readlen" \
    -p "$threads" \
    $more_args


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
echo "# Done with script"
date


# INFO -------------------------------------------------------------------------
#> RSEM-EVAL is a reference-free evaluation method based on a novel probabilistic
#> model that depends only on an assembly and the RNA-Seq reads used for its construction.
#> Unlike N50, RSEM-EVAL combines multiple factors, including the compactness of an assembly
#> and the support of the assembly from the RNA-Seq data, into a single, statistically-principled evaluation score.
#> This score can be used to select a best assembler, optimize an assembler's parameters, and guide new assembler design as an objective function.
#> In addition, for each contig within an assembly, RSEM-EVAL provides a score that
#> assesses how well that contig is supported by the RNA-Seq data and can be used to filter unnecessary contigs.
