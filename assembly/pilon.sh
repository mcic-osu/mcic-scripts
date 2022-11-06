#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=pilon
#SBATCH --output=slurm-pilon-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  RUN PILON TO POLISH A GENOME"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--genome     <file>  Genome assembly FASTA file"
    echo "  -I/--bam        <file>  BAM file of Illumina reads mapped to the assembly"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  -p/--out_prefix <str>   Prefix for output files"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -a/--more_args  <string> Quoted string with additional argument(s) to pass to Pilon"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -N/--dryrun             Dry run: don't execute commands, only parse arguments and report"
    echo "  -x/--debug              Run the script in debug mode (print all code)"
    echo "  -v/--version            Print the version of Pilon and exit"
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
    echo "    - https://github.com/broadinstitute/pilon/wiki"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/pilon-1.24
}

## Print version
Print_version() {
    Load_software
    pilon --version
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
debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
genome=""
bamdir="" && bam_arg=""
outdir=""
more_args=""


## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --genome )         shift && genome=$1 ;;
        -I | --bamdir )         shift && bamdir=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -p | --out_prefix )     shift && out_prefix=$1 ;;
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
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## Build BAM arg
for bam in "$bamdir"/*bam; do
    bam_arg="$bam_arg --frags $bam"
done

## Bash script settings
set -euo pipefail

## Check input
[[ $genome = "" ]] && Die "Please specify an input genome FASTA file with -i/--genome"
[[ $bamdir = "" ]] && Die "Please specify an input BAM dir with -I/--bam"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o/--outdir"
[[ $out_prefix = "" ]] && Die "Please specify an output prefix with -p/--out_prefix"
[[ ! -f $genome ]] && Die "Genome FASTA $genome does not exist"
[[ ! -d $bamdir ]] && Die "BAM dir $bamdir does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT PILON.SH"
date
echo "=========================================================================="
echo "Input genome assembly FASTA:  $genome"
echo "Input BAM:                    $bam"
echo "Output dir:                   $outdir"
echo "Output prefix:                $out_prefix"
[[ $more_args != "" ]] && echo "Other arguments for Pilon:    $more_args"
echo "Listing input genome FASTA:"
ls -lh "$genome"
echo "BAM file argument:"
echo "$bam_arg"
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Create the output directory
    mkdir -p "$outdir"/logs
fi

echo -e "\n## Now running Pilon..."
[[ "$debug" = false ]] && set -o xtrace

pilon \
    --genome "$genome" \
    $bam_arg \
    --outdir "$outdir" \
    --prefix "$out_prefix" \
    --fix bases \
    --threads "$n_threads"
    $more_args

[[ "$debug" = false ]] && set +o xtrace


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
