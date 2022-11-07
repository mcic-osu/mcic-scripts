#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=medaka
#SBATCH --output=slurm-medaka-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              RUN MEDAKA TO POLISH A GENOME ASSEMBLY"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--reads      <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  -r/--assembly   <file>  Input assembly: FASTA file (to be corrected)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  -m/--model      <str>   Medaka model, see the Medaka docs at https://github.com/nanoporetech/medaka#models"
    echo "                          Get a full list of possible models by running:"
    echo "                              module load miniconda3/4.12.0-py39"
    echo "                              source activate /fs/ess/PAS0471/jelmer/conda/medaka-1.7.2"
    echo "                              medaka tools list_models"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args  <str>      Quoted string with additional argument(s) to pass to Medaka"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of Medaka and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i TODO -o results/TODO "
    echo "  sbatch $0 -i TODO -o results/TODO -a \"-x TODO\""
    echo
    echo "HARDCODED PARAMETERS:"
    echo "    - ..."
    echo
    echo "OUTPUT:"
    echo "    - ..."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "    - https://github.com/nanoporetech/medaka"
    echo
}

## Load software
Load_software() {
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate; done
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/medaka-1.7.2
}

## Print version
Print_version() {
    Load_software
    medaka_consensus --version
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
## Constants

## Option defaults
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
reads=""
assembly=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )          shift && reads=$1 ;;
        -r | --assembly )       shift && assembly=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -m | --model )          shift && model=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true ;;
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
[[ $reads = "" ]] && Die "Please specify a file with input reads with -i"
[[ $assembly = "" ]] && Die "Please specify an input assembly with -r"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"
[[ ! -f $assembly ]] && Die "Input assembly file $assembly does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT MEDAKA.SH"
date
echo "=========================================================================="
echo "Input reads (FASTQ) file:         $reads"
echo "Input assembly (FASTA) file:      $assembly"
echo "Output dir:                       $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Medaka:    $more_args"
echo
echo "Listing input files:"
ls -lh "$reads" "$assembly" 
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs

    ## Run
    echo -e "\n## Now running Medaka..."
    [[ "$debug" = false ]] && set -o xtrace
    
    medaka_consensus \
        -i "$reads" \
        -d "$assembly" \
        -o "$outdir" \
        -t "$threads" \
        -m "$model"

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
