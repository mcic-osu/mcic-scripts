#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trimmomatic
#SBATCH --output=slurm-trimmomatic-%j.out

#TODO - Use adapter file by default, then test the script

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  TODO FUNCTION OF THIS SCRIPT"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
        echo "  -i STRING     Input R1 (forward reads) sequence file (name of R2 will be inferred)"
    echo "  -o STRING     Output directory (will be created if needed)"
    echo "  -i/--infile     <file>  Input file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
        echo "  -a STRING     Adapter file                       [default: 'none']"
    echo "                Possible values: 'NexteraPE-PE.fa', 'TruSeq2-PE.fa', 'TruSeq3-PE.fa'"
    echo "                Or provide the path to your own FASTA file with adapters, e.g. 'adapters.fa' from BBduk"
    echo "                With the default, no adapter trimming is done."
    echo "  -A STRING     Adapter removal parameters         [default: '2:30:10:2:True']"
    echo "  -h            Print this help message and exit"
    echo "  -p STRING     Trimming parameters for Trimmomatic"
    echo "                                                   [default: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36']"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to TODO_THIS_SOFTWARE"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for TODO_THIS_SOFTWARE and exit"
    echo "  -v/--version            Print the version of TODO_THIS_SOFTWARE and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/A1_R1.fastq.gz -o results/trimmomatic -a metadata/adapters.fa"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: http://www.usadellab.org/cms/?page=trimmomatic"
    echo "  - Paper: https://academic.oup.com/bioinformatics/article/30/15/2114/2390096"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/trimmomatic-0.39
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    trimmomatic -version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    trimmomatic -h
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
R1_in=""
outdir=""
adapter_file="" && adapter_arg=""
trim_arg=""

trim_param="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
# Same as example on https://github.com/usadellab/Trimmomatic and https://rpubs.com/ednachiang/MetaG_Pipeline
# Alternatively, example of a much stricter mode: "AVGQUAL:28 LEADING:20 TRAILING:20 MINLEN:36"

adapter_param="2:30:10:2:True"
# Same as example on https://github.com/usadellab/Trimmomatic

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1_in=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --adapter_param )   shift && adapter_param=$1 ;;
        --adapter_file )    shift && adapter_file=$1 ;;
        --trim_param )      shift && trim_param=$1 ;;
        --more_args )       shift && more_args=$1 ;;
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
# Bash script settings
set -euo pipefail

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ "$R1_in"  = "" ]] && Die "Please provide an R1 input FASTQ file with -i/--R1"
[[ "$outdir"  = "" ]] && Die "Please provide an output dir with -o/--outdir"
[[ ! -f $R1_in ]] && Die "Input file R1_in ($R1_in) does not exist"

# Process parameters
file_ext=$(basename "$R1_in" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1_in" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")
R1_basename=$(basename "$R1_in" "$file_ext")
R2_basename=$(basename "$R2_in" "$file_ext")

# Define output files
discard_dir="$outdir"/discard                          # Dir for discarded sequences
trimstats_file="$outdir"/logs/"$sample_id".trimstats.txt # File with Trimmomatic stdout
R1_out="$outdir"/"$R1_basename".fastq.gz # Output R1 file
R2_out="$outdir"/"$R2_basename".fastq.gz # Output R2 file
R1_discard=$discard_dir/"$R1_basename"_U1.fastq.gz # Output file for discarded R1 reads
R2_discard=$discard_dir/"$R2_basename"_U2.fastq.gz # Output file for discarded R2 reads

# Adapter parameters argument - As in the example here https://github.com/usadellab/Trimmomatic
[[ $adapter_file != "" ]] && adapter_arg=" ILLUMINACLIP:$adapter_file:$adapter_param"
[[ $trim_param != "" ]] && trim_arg=" $trim_param"

# Check parameters
[[ ! -f $R2_in ]] && Die "Input file R1_in ($R2_in) does not exist"
[[ "$R1_in" = "$R2_in" ]] && echo "Input R1 and R2 FASTQ files are the same file: $R1"
[[ "$R1_in" = "$R1_out" ]] && echo "Input R1 and output R1 FASTQ files are the same file: $R1"
[[ "$R2_in" = "$R2_out" ]] && echo "Input R2 and output R2 FASTQ files are the same file: $R1"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TRIMMOMATIC.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input R1 file:                    $R1_in"
echo "Input R2 file:                    $R2_in"
echo "Output dir:                       $outdir"
echo "Output R1 file:                   $R1_out"
echo "Output R2 file:                   $R2_out"
[[ "$trim_arg" != "" ]] && echo "Trimming argument:                $trim_arg"
[[ "$adapter_arg" != "" ]] && echo "Adapter argument:                 $adapter_arg"
[[ $more_args != "" ]] && echo "Other arguments for Trimmomatic:$more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$R1_in" "$R2_in"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create output dirs
echo -e "\n# Creating the output directories..."
${e}mkdir -p "$outdir/logs" "$discard_dir"

# Run Trimmomatic
echo -e "\n# Starting the Trimmomatic run..."
${e}Time \
    trimmomatic PE \
    -threads "$threads" \
    "$R1_in" "$R2_in" \
    "$R1_out" "$R1_discard" \
    "$R2_out" "$R2_discard"${adapter_arg}${trim_arg}${more_args} \
    2>&1 | tee "$trimstats_file"

# Count the number of reads
echo -e "\n# Now counting the number of reads in the in- and output..."
nreads_raw=$(zcat "$R1_in" | awk '{s++}END{print s/4}')
nreads_trim=$(zcat "$R1_out" | awk '{s++} END{print s/4}')
echo -e "Number of raw / trimmed read-pairs: $nreads_raw / $nreads_trim"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lh "$PWD"/"$outdir"/"$R1_out" "$PWD"/"$outdir"/"$R2_out"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
echo
