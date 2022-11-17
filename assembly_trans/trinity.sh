#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=170G
#SBATCH --cpus-per-task=42
#SBATCH --job-name=trinity
#SBATCH --output=slurm-trinity-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "============================================================================"
    echo "                            $0"
    echo "                     Run Trinity to assemble a transcriptome"
    echo "============================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i                          <file>  Input dir with FASTQ files (de novo assembly) OR a BAM file (genome-guided assembly)"
    echo "  -o/--outdir                 <dir>   Output dir (will be created if needed)"
    echo "                                      NOTE: The output directory needs to include 'trinity' in its name"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --SS_lib_type / --strandedness <str>   RNAseq library type: 'RF'/'reverse', 'FR'/'forward', or 'unstranded'   [default: 'RF']"
    echo "  --min_contig_length         <int>   Minimum contig length                   [default: 200 (= Trinity default)]"
    echo "  --normalize                 <bool>  Whether to normalize reads, 'true' or 'false' [default: 'true']"
    echo "  --normalize_max_read_cov    <int>   Normalize to this coverage              [default: 200 (= Trinity default)]"
    echo "  --genome_guided                     Genome-guided assembly                  [default: de novo assembly]"
    echo "  --genome_guided_max_intron  <int>   Max intron size for genome-guided assembly  [default: 10000]"
    echo "  --more-args                 <str>   Quoted string with additional argument(s) to pass to Trinity"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                            Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                             Run the script in debug mode (print all code)"
    echo "  -h                                  Print this help message and exit"
    echo "  --help                              Print the help for Trinity and exit"
    echo "  -v/--version                        Print the version of Trinity and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/ -o results/trinity"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/trinityrnaseq/trinityrnaseq/wiki"
    echo
}

## Load the software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/trinity-2.13.2
}

## Print version
Print_version() {
    Load_software
    Trinity --version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    Trinity --help
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
#                     CONSTANTS AND DEFAULTS
# ==============================================================================
strandedness="RF"
genome_guided=false
genome_guided_max_intron=10000  # Only for genome-guided assembly
normalize=true && normalize_arg=""
min_contig_length=200
normalize_max_read_cov=200

debug=false
dryrun=false && e=""
slurm=true

# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
## Placeholder defaults
indir=""
outdir=""
bam=""
more_args=""
mem_gb=4   # Will be changed if this is a SLURm job

## Parse command-line options
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i )                            shift && input=$1 ;;
        -o | --outdir )                 shift && outdir=$1 ;;
        --SS_lib_type | --strandedness) shift && strandedness=$1 ;;
        --min_contig_length )           shift && min_contig_length=$1 ;;
        --normalize )                   shift && normalize=$1 ;;
        --normalize_max_read_cov )      shift && normalize_max_read_cov=$1 ;;
        --genome_guided )               genome_guided=true ;;
        --genome_guided_max_intron )    shift && genome_guided_max_intron=$1 ;;
        --more-args )                   shift && more_args=$1 ;;
        -v | --version )                Print_version; exit 0 ;;
        -h )                            Print_help; exit 0 ;;
        --help )                        Print_help_program; exit 0;;
        --dryrun )                      dryrun=true && e="echo ";;
        --debug )                       debug=true ;;
        * )                             Die "Invalid option $1" "$all_args" ;;
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

## Check input
[[ "$input" = "" ]] && Die "Please specify input with -i" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ "$normalize" != "true" && "$normalize" != "false" ]] && Die "--normalize should be 'true' or 'false', instead it is $normalize"

## Create comma-delimited list of FASTQ files:
if [[ "$genome_guided" = false ]]; then
    indir="$input"
    [[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"
    R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
    R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')
else
    bam="$input"
    [[ ! -f "$bam" ]] && Die "Input BAM file $bam does not exist"
fi

## Library type argument
if [[ "$strandedness" = "unstranded" ]]; then
    strand_arg=""
elif [[ "$strandedness" = "reverse" ]]; then
    strand_arg="--SS_lib_type RF"
elif [[ "$strandedness" = "forward" ]]; then
    strand_arg="--SS_lib_type FR"
else
    strand_arg="--SS_lib_type $strandedness"
fi

## Normalization arg
[[ "$normalize" = "false" ]] && normalize_arg="--no_normalize_reads"

## Define memory in GB
[[ "$slurm" = true ]] && mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TRINITY.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input dir:                        $indir"
echo "Output dir:                       $outdir"
echo "Strandedness / strand argument:   $strandedness / $strand_arg"
echo "Minimum contig length:            $min_contig_length"
echo "Normalize reads:                  $normalize"
echo "Max cov. to normalize reads to:   $normalize_max_read_cov"
[[ $more_args != "" ]] && echo "Other arguments for Trinity:      $more_args"
echo "Number of threads/cores:          $threads"
echo "Memory in GB:                     $mem_gb"
echo
if [[ "$genome_guided" = false ]]; then
    echo "List of R1 FASTQ files:           $R1_list"
    echo "List of R2 FASTQ files:           $R2_list"
    echo
    echo "Listing the input file(s):"
    ls -lh "$indir"
else
    echo "BAM file:                         $bam"
    echo "Max intron size:                  $genome_guided_max_intron"
    echo "Listing the input file(s):"
    ls -lh "$bam"
fi

[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

# MAIN -------------------------------------------------------------------------
echo "## Starting Trinity run..."

[[ "$dryrun" = false ]] && set -o xtrace

if [[ "$genome_guided" = false ]]; then
    ${e}Time Trinity \
        --seqType fq \
        --left "$R1_list" \
        --right "$R2_list" \
        --min_contig_length "$min_contig_length" \
        --normalize_max_read_cov "$normalize_max_read_cov" \
        $normalize_arg \
        $strand_arg \
        --output "$outdir" \
        --max_memory "$mem_gb" \
        --CPU "$threads" \
        --verbose \
        $more_args
else
    ${e}Time Trinity \
        --genome_guided_bam "$bam" \
        --genome_guided_max_intron "$genome_guided_max_intron" \
        --min_contig_length "$min_contig_length" \
        --normalize_max_read_cov "$normalize_max_read_cov" \
        $normalize_arg \
        $strand_arg \
        --output "$outdir" \
        --max_memory "$mem_gb" \
        --CPU "$threads" \
        --verbose
fi

[[ "$debug" = false ]] && set +o xtrace

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
