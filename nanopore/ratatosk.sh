#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --job-name=ratatosk
#SBATCH --output=slurm-ratatosk-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         RUN RATATOSK TO CORRECT LONG READS WITH ILLUMINA READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--fq_long       <file>     Input long-read FASTQ file"
    echo "    -I/--fq_short_list <file>     List with input short-read FASTQ file(s)"
    echo "    -o/--outdir        <dir>      Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    --insert_size      <integer>  Insert size - set to read length for single-end reads  [default: 500]"
    echo "    -a/--more_args     <string>   Quoted string with additional argument(s) to pass to Ratatosk"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h/--help                     Print this help message and exit"
    echo "    -N/--dryrun                   Dry run: don't execute commands, only parse arguments and report"
    echo "    -x/--debug                    Run the script in debug mode (print all code)"
    echo "    -v/--version                  Print the version of Ratatosk and exit"
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
#Load_software() {
    #module load miniconda3/4.12.0-py39
    #source activate /fs/ess/PAS0471/jelmer/conda/fmlrc2-0.1.7
    #! Currently installed in /users/PAS0471/jelmer/.local/bin/
#}

## Print version
Print_version() {
    #Load_software
    Ratatosk --version
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
insert_size=500

debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
fq_long=""
fq_short_list=""
outdir=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --fq_long )        shift && fq_long=$1 ;;
        -I | --fq_short_list )  shift && fq_short_list=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --insert_size )         shift && insert_size=$1 ;;
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
#[[ "$dryrun" = false ]] && Load_software

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
[[ $fq_long = "" ]] && Die "Please specify a long-read FASTQ file with --fq_long"
[[ $fq_short_list = "" ]] && Die "Please specify a list with short-read FASTQ files with --fq_short_list"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $fq_long ]] && Die "Input file $fq_long does not exist"
[[ ! -f $fq_short_list ]] && Die "Input file $fq_short_list does not exist"

## Define output files (NOTE: Ratatosk will add .fastq)
fq_out="$outdir"/$(basename "$fq_long" .fastq.gz)

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RATATOSK.SH"
date
echo "=========================================================================="
echo "Input long-read FASTQ file:               $fq_long"
echo "List with input short-read FASTQ file(s): $fq_short_list"
echo "Output dir:                               $outdir"
echo "Insert size:                              $insert_size"
echo
echo "Output FASTQ file:                        $fq_out"
[[ $more_args != "" ]] && echo "Other arguments for Ratatosk:    $more_args"
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
    echo -e "\n## Now running Ratatosk..."
    [[ "$debug" = false ]] && set -o xtrace
    Ratatosk correct \
        -v \
        --cores "$threads" \
        --insert-sz "$insert_size" \
        --in-short "$fq_short_list" \
        --in-long "$fq_long" \
        --out-long "$fq_out"
    [[ "$debug" = false ]] && set +o xtrace

    ## Gzip FASTQ file
    echo -e "\n## Now gzipping the output FASTQ file..."
    gzip "$fq_out".fastq
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
