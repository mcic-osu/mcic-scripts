#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=43
#SBATCH --mem=172G
#SBATCH --job-name=flye
#SBATCH --output=slurm-flye-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "           RUN FLYE TO ASSEMBLE A GENOME WITH LONG READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--infile FILE       Input file"
    echo "    -o/--outdir DIR        Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -g/--genome_size STRING  Genome size estimate, e.g '4.6m'         [default: no estimate]"
    echo "    -a/--more_args STRING  Quoted string with additional argument(s) to pass to Flye"
    echo "    -r/--resume            Resume previous run"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h/--help              Print this help message and exit"
    echo "    -N/--dryrun            Dry run: don't execute commands, only parse arguments and report"
    echo "    -x/--debug             Run the script in debug mode (print all code)"
    echo "    -v/--version           Print the version of Flye and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "    sbatch $0 -i TODO -o results/TODO "
    echo "    sbatch $0 -i TODO -o results/TODO -a \"-x TODO\""
    echo
    echo "HARDCODED PARAMETERS:"
    echo "    - ..."
    echo
    echo "OUTPUT:"
    echo "    - ..."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - ..."
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/flye-2.9.1
}

## Print version
Print_version() {
    Load_software
    flye --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
resume=false && resume_arg=""
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
outdir=""
genome_size="" && genome_size_arg=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -s | --genome_size )    shift && genome_size=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        -r | --resume )         resume=true ;;
        -X | --debug )          debug=true ;;
        -N | --dryrun )         dryrun=true ;;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_help; Die "Invalid option $1" ;;
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
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## Build other args
[[ "$genome_size" != "" ]] && genome_size_arg="--genome-size $genome_size"
[[ "$resume" = true ]] && resume_arg="--resume"

## Bash script settings
set -euo pipefail

## Check input
[[ $infile = "" ]] && Die "Please specify an input file with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FLYE.SH"
date
echo "=========================================================================="
echo "Input file:                  $infile"
echo "Output dir:                  $outdir"
echo "Genome size:                 $genome_size"
echo "Resume previous run:         $resume"
[[ $more_args != "" ]] && echo "Other arguments for Flye:    $more_args"
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "Listing input file:"
ls -lh "$infile"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs

    echo -e "\n## Now running Flye..."
    [[ "$debug" = false ]] && set -o xtrace

    flye \
        --threads "$n_threads" \
        --nano-raw "$infile" \
        --out-dir "$outdir" \
        $resume_arg \
        $genome_size_arg \
        $more_args

    [[ "$debug" = false ]] && set +o xtrace
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
