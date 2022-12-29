#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=quickmerge
#SBATCH --output=slurm-quickmerge-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              MERGE LONG-READ ASSEMBLIES WITH QUICKMERGE"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -o <output file> --query <input-assembly1> --ref <input-assembly2> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--merged     <file>  Merged assembly FASTA file"
    echo "  --query         <file>  Input assembly FASTA #1, used as 'query'"
    echo "                          Note: Quickmerge will use --ref to improve the --query, so the output will be most like the --query"
    echo "                          Therefore, use the best assembly as the --query (see https://github.com/mahulchak/quickmerge/wiki)."
    echo "  --ref           <file>  Input assembly FASTA #2, used as 'reference'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --minlen_merge  <int>"
    echo "  --minlen_anchor <int>"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Quickmerge"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Quickmerge and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -o results/quickmerge.fasta --query results/assembly1.fasta --ref results/assembly2.fasta"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/mahulchak/quickmerge"
    echo "  - More docs: https://github.com/mahulchak/quickmerge/wiki"
    echo "  - Paper: https://academic.oup.com/nar/article/44/19/e147/2468393"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/quickmerge-env
    set -u
}

# Print help for the focal program
Print_help_program() {
    Load_software
    merge_wrapper.py --help
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
minlen_merge=10000     # Quickmerge default is 5000 but it says higher values are recommended
minlen_anchor=100000   #TODO auto-determine this value
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
query=""
ref=""
merged=""
minlen_anchor=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --merged )     shift && merged=$1 ;;
        --query )           shift && query=$1 ;;
        --ref )             shift && ref=$1 ;;
        --minlen_merge )    shift && minlen_merge=$1 ;;
        --minlen_anchor )   shift && minlen_anchor=$1 ;;
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

# Load software
[[ "$dryrun" = false ]] && Load_software

# Check input
[[ "$query" = "" ]] && Die "Please specify an input 'query' file with --query" "$all_args"
[[ "$ref" = "" ]] && Die "Please specify an input 'reference' file with --ref" "$all_args"
[[ "$merged" = "" ]] && Die "Please specify an output assembly with -o/--merged" "$all_args"
[[ ! -f "$query" ]] && Die "Input query file $query does not exist"
[[ ! -f "$ref" ]] && Die "Input reference file $ref does not exist"

# Make paths absolute (we have to move into the outdir)
[[ ! "$query" =~ ^/ ]] && query="$PWD"/"$query"
[[ ! "$ref" =~ ^/ ]] && ref="$PWD"/"$ref"
[[ ! "$merged" =~ ^/ ]] && merged="$PWD"/"$merged"

# Determine the output dir
outdir=$(dirname "$merged")

# Report
echo
echo "=========================================================================="
echo "                 STARTING SCRIPT QUICKMERGE.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Query assembly (input 1):             $query"
echo "Reference assembly (input 2):         $ref"
echo "Merged assembly (output):             $merged"
echo "Min. anchor length (Quickmerge -l):   $minlen_anchor"
echo "Min. merging length (Quickmerge -ml): $minlen_merge"
[[ $more_args != "" ]] && echo "Other arguments for Quickmerge:       $more_args"
echo
echo "Listing the input file(s):"
ls -lh "$query" "$ref"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
${e}mkdir -pv "$outdir"/logs

# Move into outdir (Quickmerge will create files in working dir)
cd "$outdir" || exit 

# Run
echo -e "\n# Now running Quickmerge..."
${e}Time \
    merge_wrapper.py \
        -ml "$minlen_merge" \
        -l "$minlen_anchor" \
        $more_args \
        "$query" \
        "$ref" \
        > "$merged"

#? -l: controls the length cutoff for anchor contigs.
#? A good rule of thumb is to start with the N50 of the self assembly.
#? E.g. if the N50 of your self assembly is 2Mb then use 2000000 as your cutoff.
#? Lowering this value may lead to more merging but may increase the probability of mis-joins.

#? -ml: controls the minimum alignment length to be considered for merging.
#? This is especially helpful for repeat-rich genomes. Default is 5000 but higher values (>5000) are recommended.


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
