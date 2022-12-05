#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=samtools
#SBATCH --output=slurm-samtools-%j.out

# ==============================================================================
#                           FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                       RUN SAMTOOLS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> --command <samtools command>"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <dir>       Input file"
    echo "  --command       <string>    Quoted string with additional argument(s) to pass to Samtools"
    echo
    echo "UTILITY OPTIONS:" 
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v/--version                Print the version of Samtools and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  # To extract unmapped reads to a FASTQ file:"
    echo "  sbatch $0 -i results/my.bam --command 'fastq -n -f 4 -0 results/my.fastq.gz'"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - http://www.htslib.org/doc/samtools.html"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/samtools
}

## Print version
Print_version() {
    Load_software
    samtools --version | head -n 2
}

## Print args
Print_args() {
    echo -e "\n# Arguments passed to the script:"
    echo "$*"
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo "For help, run this script with the '-h' / '--help' option"
    echo -e "Exiting\n" >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
command=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        --command )             shift && command=$1 ;;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        --dryrun )              dryrun=true ;;
        --debug )               debug=true ;;
        * )                     Print_args "$all_args"; Die "Invalid option $1" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## Get number of threads
if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
    threads="$SLURM_CPUS_PER_TASK"
elif [[ -n "$SLURM_NTASKS" ]]; then
    threads="$SLURM_NTASKS"
else
    threads=1
fi

## Bash script settings
set -euo pipefail

## Check input
[[ $infile = "" ]] && Print_args "$all_args" && Die "Please specify an input file with -i"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT SAMTOOLS.SH"
date
echo "=========================================================================="
echo "Input dir:                    $infile"
echo "Command:                      $command"
echo "Listing the input file:"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN\n"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    
    ## Run
    echo -e "\n## Now running Samtools..."
    [[ "$debug" = false ]] && set -o xtrace\
    
    samtools \
        $command \
        --threads $threads \
        "$infile"

    [[ "$debug" = false ]] && set +o xtrace

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    echo
    echo "========================================================================="
    echo "## Version used:"
    Print_version
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
