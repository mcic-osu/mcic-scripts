#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
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
    echo "  sbatch $0 -i <input FASTQ> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile   <file>    Input FASTQ file"
    echo "  -o/--outdir   <dir>     Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --genome_size <str>     Genome size estimate, e.g '4.6m'            [default: no estimate]"
    echo "  --iterations  <int>     Number of polishing iterations              [default: 1]"
    echo "  --more_args   <str>     Quoted string with additional argument(s) to pass to Flye"
    echo "  --resume                Resume previous run"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of Flye and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/minion/my.fastq -o results/flye"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - ..."
    echo
    echo "OUTPUT:"
    echo "  - ..."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md"
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
iterations=1
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
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --genome_size )     shift && genome_size=$1 ;;
        --iterations )      shift && iterations=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        --resume )          resume=true ;;
        -v | --version )    Print_version; exit ;;
        -h | --help )       Print_help; exit ;;
        --debug )           debug=true ;;
        --dryrun )          dryrun=true ;;
        * )                 Print_help; Die "Invalid option $1" ;;
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
echo "Nr of polishing iterations:  $iterations"
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
        --iterations "$iterations" \
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
    ls -lhd "$PWD"/"$outdir"/*
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
