#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=bactopia_datasets
#SBATCH --output=slurm-bactopia-datasets-%j.out


# PARSE COMMAND-LINE OPTIONS ---------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "$0: Download Bactopia datasets/databases"
    echo "=================================================================================================="
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -o DIR      Output directory (will be created if needed)"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -s STRING   Species name"
    echo "    -a STRING   Additional arguments to pass to 'bactopia datasets'"
    echo "    -h          Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "    sbatch $0 -i data/meta/samplesheet.tsv -o results/bactopia_datasets"
    echo
    echo "OUTPUT:"
    echo "------------------"
    echo "      - ...."
    echo
    echo "BACTOPIA DOCUMENTATION:"
    echo "------------------"
    echo " - https://bactopia.github.io/ "
}

## Option defaults
species=""
more_args=""
debug=false

## Parse command-line options
while getopts 'o:s:a:hx' flag; do
    case "${flag}" in
        o) outdir="$OPTARG" ;;
        s) species="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        x) debug=true ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SOFTWARE SETUP ---------------------------------------------------------------
## Load Conda environment
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bactopia


# OTHER SETUP ------------------------------------------------------------------
## Bash strict settings
set -ueo pipefail

## Check input
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" && exit 1

## Get nr of threads
set +u
if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
    n_threads="$SLURM_NTASKS"
else
    n_threads="$SLURM_CPUS_PER_TASK"
fi
set -u

## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT BACTOPIA-DATASETS.SH"
date
echo -e "==========================================================================\n"
echo "## Output dir:                        $outdir"
echo "## Species:                           $species"
[[ "$more_args" != "" ]] && echo "## Additional arguments:              $more_args"
echo -e "-------------------------\n"


# MAIN -------------------------------------------------------------------------
## Make necessary dirs
mkdir -p "$outdir"

## Define the workflow command
command="bactopia datasets \
            --outdir $outdir \
            --species '$species' \
            --include_genus \
            --cpus $n_threads \
            $more_args"
echo "## Running 'bactopia datasets' using the following command:"
echo
echo "$command"
echo
[[ "$debug" = false ]] && eval "$command"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script bactopia-datasets.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
