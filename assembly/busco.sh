#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --job-name=busco
#SBATCH --output=slurm-busco-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         Run BUSCO to check a transcriptome or genome assembly"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-FASTA> -o <output-dir> -d <db-name> ..."
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>      Input assembly FASTA file"
    echo "  -o/--outdir     <dir>       Output dir (will be created if needed)"
    echo "  -d/--db         <string>    Busco database name (see https://busco.ezlab.org/list_of_lineages.html)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -m/--mode       <string>    Mode, i.e. assembly type               [default: 'genome']"
    echo "                              Valid options: 'genome', 'transcripttome', or 'proteins'"
    echo "  --more-args     <str>       Quoted string with additional argument(s) to pass to Busco"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Busco and exit"
    echo "  -v/--version            Print the version of Busco and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly.fa -o results/busco -d bacteria_odb"
    echo
    echo "DOCUMENTATION:"
    echo "  - https://busco.ezlab.org/busco_userguide.html"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/busco
}

## Print version
Print_version() {
    Load_software
    busco --version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    busco --help
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
assembly_type=genome

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
outdir=""
busco_db=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        -d | --db )         shift && busco_db=$1 ;;
        -m | --mode )       shift && assembly_type=$1 ;;
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

## If needed, make input path absolute because we have to move into the outdir
[[ ! $infile =~ ^/ ]] && infile="$PWD"/"$infile"

## Get a sample/assembly ID from the filename
fileID=$(basename "$infile" | sed -E 's/.fn?as?t?a?//')

## Check input
[[ "$infile" = "" ]] && Die "Please specify an input FASTA file with -i"
[[ "$outdir" = "" ]] && Die "Please specify an output dir file with -o"
[[ "$busco_db" = "" ]] && Die "Please specify a Busco database name with -d"
[[ ! -f "$infile" ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT BUSCO.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input FASTA:                      $infile"
echo "Output dir:                       $outdir"
echo "BUSCO db:                         $busco_db"
echo "Mode (assembly type):             $assembly_type"
[[ $more_args != "" ]] && echo "Other arguments for Busco:        $more_args"
echo
echo "Number of cores:                  $threads"
echo "Assembly ID (inferred):           $fileID"
echo "Listing the input file(s):"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create output dir if needed
    mkdir -p "$outdir"/logs

    ## Move into output dir
    cd "$outdir" || exit 1

    ## Run BUSCO
    echo "## Now running BUSCO..."
    [[ "$debug" = false ]] && set -o xtrace

    ${e}Time busco \
        -i "$infile" \
        -o "$fileID" \
        -l "$busco_db" \
        -m "$assembly_type" \
        -c "$threads" \
        --force
    
    [[ "$debug" = false ]] && set +o xtrace

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
