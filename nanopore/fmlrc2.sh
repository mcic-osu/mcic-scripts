#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --job-name=fmlrc2
#SBATCH --output=slurm-fmlrc2-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         RUN FMLRC2 TO CORRECT LONG READS WITH ILLUMINA READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--fq_long    <file>    Input long-read FASTQ file"
    echo "    -I/--fq_short   <string>  Input short-read FASTQ file(s)"
    echo "    -o/--outdir     <dir>     Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -a/--more_args  <string>  Quoted string with additional argument(s) to pass to FMLRC2"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h/--help                 Print this help message and exit"
    echo "    -N/--dryrun               Dry run: don't execute commands, only parse arguments and report"
    echo "    -x/--debug                Run the script in debug mode (print all code)"
    echo "    -v/--version              Print the version of FMLRC2 and exit"
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
    source activate /fs/ess/PAS0471/jelmer/conda/fmlrc2-0.1.7
}

## Print version
Print_version() {
    Load_software
    fmlrc2 --version
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
fq_long=""
fq_short=""
outdir=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --fq_long )        shift && fq_long=$1 ;;
        -I | --fq_short )       shift && fq_short=$1 ;;
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
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## Bash script settings
set -euo pipefail

## Check input
[[ $fq_long = "" ]] && Die "Please specify a long-read FASTQ file with --fq_long"
[[ $fq_short = "" ]] && Die "Please specify one more short-read FASTQ files with --fq_short"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $fq_long ]] && Die "Input file $fq_long does not exist"

## Define output files
bwt_out="$outdir"/comp_msbwt.npy
fasta_out="$outdir"/$(basename "$fq_long" .fastq.gz).fasta

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FMLRC2.SH"
date
echo "=========================================================================="
echo "Input long-read FASTQ file:     $fq_long"
echo "Input short-read FASTQ file(s): $fq_short"
echo "Output dir:                     $outdir"
echo
echo "Output BWT file:                $bwt_out"
echo "Output FASTA file:              $fasta_out"
[[ $more_args != "" ]] && echo "Other arguments for FMLRC2:    $more_args"
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

[[ "$debug" = false ]] && set -o xtrace

echo -e "\n## Now building the BWT with ropebwt2..."
gunzip -c $fq_short | \
    awk 'NR % 4 == 2' | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert "$bwt_out"

echo -e "\n## Now running fmlrc2..."
fmlrc2 \
    --threads "$n_threads" \
    --cache_size 11 \
    $more_args \
    "$bwt_out" \
    "$fq_long" \
    "$fasta_out"

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
