#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=racon
#SBATCH --output=slurm-racon-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              RUN RACON TO POLISH A GENOME ASSEMBLY"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--reads      <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  -b/--align      <file>  Input aligments: BAM/PAF file (reads mapped to assembly)"
    echo "  -r/--assembly   <file>  Input assembly: FASTA file (to be corrected)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -a/--more_args  <string> Quoted string with additional argument(s) to pass to Racon"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -N/--dryrun             Dry run: don't execute commands, only parse arguments and report"
    echo "  -x/--debug              Run the script in debug mode (print all code)"
    echo "  -v/--version            Print the version of Racon and exit"
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
    echo "    - ..."
    echo
}

## Load software
Load_software() {
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate; done
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/racon-1.5.0
}

## Print version
Print_version() {
    Load_software
    racon --version
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
align=""
assembly=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )          shift && reads=$1 ;;
        -b | --align )          shift && align=$1 ;;
        -r | --assembly )       shift && assembly=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
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
if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
    threads="$SLURM_CPUS_PER_TASK"
elif [[ -n "$SLURM_NTASKS" ]]; then
    threads="$SLURM_NTASKS"
else
    threads=1
fi

## Bash script settings
set -euo pipefail

## Define output file
assembly_ext=$(echo "$assembly" | sed -E 's/.*(\.fn?a?s?t?a$)/\1/')
outfile="$outdir"/$(basename "$assembly" "$assembly_ext").fasta

## Check input
[[ $reads = "" ]] && Die "Please specify a file with input reads with -i"
[[ $align = "" ]] && Die "Please specify a file with input alignments with -a"
[[ $assembly = "" ]] && Die "Please specify an input assembly with -r"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"
[[ ! -f $align ]] && Die "Input align file $align does not exist"
[[ ! -f $assembly ]] && Die "Input assembly file $assembly does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RACON.SH"
date
echo "=========================================================================="
echo "Input reads (FASTQ) file:         $reads"
echo "Input alignment (BAM) file:       $align"
echo "Input assembly (FASTA) file:      $assembly"
echo "Output dir:                       $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Racon:    $more_args"
echo "Listing input files:"
ls -lh "$reads" "$align" "$assembly" 
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs

    ## Run
    echo -e "\n## Now running Racon..."
    [[ "$debug" = false ]] && set -o xtrace
    
    racon \
        "$reads" \
        "$align" \
        "$assembly" \
        --threads "$threads"
        $more_args \
        > "$outfile"

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
