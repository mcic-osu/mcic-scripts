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
    echo "  -i STRING         Input FASTA file with the assembly"
    echo "  -o STRING         Output directory"
    echo "  -d STRING         Busco database name (see https://busco.ezlab.org/list_of_lineages.html)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -m STRING         Run mode, i.e. assembly type               [default: 'genome']"
    echo "                    Valid options: 'genome', 'transcripttome', or 'proteins'"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly/assembly.fa -o results/BUSCO -d bacteria_odb"
    echo
    echo "DOCUMENTATION:"
    echo "  - https://busco.ezlab.org/busco_userguide.html"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /users/PAS0471/jelmer/miniconda3/envs/busco-env
}

## Print version
Print_version() {
    Load_software
    busco --version
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
assembly_type=genome

debug=false
dryrun=false

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
fa_in=""
outdir=""
busco_db=""

## Parse command-line args
while getopts ':i:o:d:m:h' flag; do
    case "${flag}" in
        i) fa_in="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) busco_db="$OPTARG" ;;
        m) assembly_type="$OPTARG" ;;
        h) Print_help; exit 0 ;;
        \?) Die "Invalid option -$OPTARG" ;;
        :) Die "Option -$OPTARG requires an argument." ;;
    esac
done

## Check input
[[ "$fa_in" = "" ]] && Die "Please specify an input FASTA file with -i"
[[ "$outdir" = "" ]] && Die "Please specify an output dir file with -o"
[[ "$busco_db" = "" ]] && Die "Please specify a Busco database name with -d"
[[ ! -f "$fa_in" ]] && Die "Input file $fa_in does not exist"


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Software
Load_software

## Bash strict mode
set -euo pipefail

## Get number of threads
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## If needed, make input path absolute because we have to move into the outdir
[[ ! $fa_in =~ ^/ ]] && fa_in="$PWD"/"$fa_in"

## Get a sample/assembly ID from the filename
fileID=$(basename "$fa_in" | sed -E 's/.fn?as?t?a?//')

## Create output dir if needed
mkdir -p "$outdir"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT BUSCO.SH"
date
echo "=========================================================================="
echo "Input FASTA:               $fa_in"
echo "Output dir:                $outdir"
echo "BUSCO db:                  $busco_db"
echo "Mode (assembly type):      $assembly_type"
echo
echo "Number of cores:           $n_threads"
echo "Assembly ID (inferred):    $fileID"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
## Move into output dir
cd "$outdir" || exit 1

echo "## Now running busco..."

[[ "$debug" = false ]] && set -o xtrace
busco \
    -i "$fa_in" \
    -o "$fileID" \
    -l "$busco_db" \
    -m "$assembly_type" \
    -c "$n_threads" \
    --force
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
