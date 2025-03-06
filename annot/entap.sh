#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=entap
#SBATCH --output=slurm-entap-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "     Run EnTAP version 0.10.8 to annotate a transcriptome assembly"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -c <config file> -o <output dir> [ --db_dir <dir> | <db-DIAMOND-1> [<db-DIAMOND-2> ...]]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly   <file>  Input assembly (nucleotide) FASTA file"
    echo "  -c/--config     <file>  Input EnTAP config file"
    echo "                            This file can be generated by running 'entap_config.sh' (present in same dir as this script)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "DIAMOND REFERENCE DATABASES (ALSO REQUIRED):"
    echo "  To specify DIAMOND database files, specify --db_dir OR pass files as positional arguments."
    echo "  You can use the following pre-built databases:"
    echo "    - RefSeq Complete:    /fs/project/PAS0471/jelmer/refdata/entap/bin/refseq_complete.dmnd"
    echo "    - UniProt:            /fs/project/PAS0471/jelmer/refdata/entap/bin/uniprot_sprot.dmnd"
    echo "  --db_dir        <dir>   Directory with DIAMOND database files ('.dmnd')"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --bam           <file>  Input BAM filename"
    echo "                            Default is to run without a BAM file and therefore without expression level filtering"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to EnTAP"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for EnTAP and exit"
    echo "  -v/--version            Print the version of EnTAP and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -c entap_config.ini -d uniprot_sprot.fa -o results/entap"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - EnTAP documentation: https://entap.readthedocs.io/en/latest/"
    echo "  - EnTAP paper: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13106"
    echo
}

# Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/entap-0.10.8
}

# Print version
Print_version() {
    set +u
    Load_software
    EnTAP --version
    set -u
}

# Print help for the focal program
Print_help_program() {
    Load_software
    EnTAP --help
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
    echo "Memory (per node):    $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
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
# Other option defaults
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly=""
declare -a dbs
db_dir=""
db_arg=""
outdir=""
config=""
bam="" && align_arg=""
more_args=""

# Parse command-line args
all_args="$*"
count=0
while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly )   shift && assembly=$1 ;;
        -c | --config )     shift && config=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --db-dir )          shift && db_dir=$1 ;;
        --bam )             shift && bam=$1 ;;
        --more-args )       shift && more_args=$1 ;;
        -v | -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 dbs[count]=$1 && count=$(( count + 1 )) ;;
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

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Bash script settings
set -euo pipefail

# Check input
[[ "$assembly" = "" ]] && Die "Please specify an assembly file with -i/--assembly" "$all_args"
[[ "$config" = "" ]] && Die "Please specify a config file with -c/--config" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$assembly" ]] && Die "Assembly file $assembly does not exist" "$all_args"
[[ ! -f "$config" ]] && Die "Config file $config does not exist" "$all_args"
[[ "$bam" != "" ]] && [[ ! -f "$bam" ]] && Die "BAM file $bam does not exist" "$all_args"

# Make paths absolute, or EnTap will fail
[[ ! $config =~ ^/ ]] && config="$PWD"/"$config"
[[ ! $assembly =~ ^/ ]] && assembly="$PWD"/"$assembly"
[[ ! $outdir =~ ^/ ]] && outdir="$PWD"/"$outdir"

# BAM arg
if [[ "$bam" != "" ]]; then
    bam="$PWD"/"$bam"  # Make path absolute
    align_arg="--align $bam"
fi

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT ENTAP.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Input assembly:                       $assembly"
echo "Output dir:                           $outdir"
echo "Config file:                          $config"
echo "Number of entries in the assembly:    $(grep -c "^>" "$assembly")"
echo "Number of threads/cores:              $threads"
[[ "$bam" != "" ]] && echo "Input BAM file:                       $bam"
[[ $more_args != "" ]] && echo "Other arguments for EnTAP:            $more_args"

# Build database argument for EnTap
echo "# Listing the DIAMOND databases:"
if [[ "$db_dir" != "" ]]; then
    if compgen -G "$db_dir"/*dmnd > /dev/null; then
        for db in "$db_dir"/*dmnd; do
            [[ ! $db =~ ^/ ]] && db="$PWD"/"$db"
            db_arg="$db_arg -d $db"
            ls -lh "$db"
        done
    else
        echo "No DIAMOND databases found in $db_dir"
    fi
elif [[ ${#dbs[@]} -gt 0 ]]; then
    db_arg="-d ${dbs// / -d }"
    for db in "${dbs[@]}"; do ls -lh "$db"; done
else
    Die "Please specify input DIAMOND databases with --dbs, --db-dir, or as positional args"
fi

# Report, part 2
echo "Input database arg:                   $db_arg"
echo
echo "# Printing the contents of the config file:"
cat -n "$config"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create output directory
echo -e "\nCreating the output directories..."
mkdir -pv "$outdir"/logs

# Copy the entap_graphing.py script to the working dir
# Entap will call it using 'python entap_graphing.py' which won't work unless it's in the current working dir
echo -e "\nCopying the graphing script to the working dir..."
graphing_script=$(whereis entap_graphing.py | sed 's/.* //')
cp -v "$graphing_script" .

echo -e "\n# Now running EnTAP..."
${e}Time EnTAP \
    --runP \
    --ini "$config" \
    -i "$assembly" \
    --out-dir "$outdir" \
    $db_arg \
    $align_arg \
    -t "$threads" \
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
    ls -lhd "$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
