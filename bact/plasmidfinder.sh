#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=plasmidfinder
#SBATCH --output=slurm-plasmidfinder-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Run PlasmidFinder to detect plasmids in a genome."
    echo
    echo "Syntax: $0 -i <genome-FASTA> -o <output-dir> ..."
    echo
    echo "Required options:"
    echo "-i FILE           Genome (nucleotide) FASTA file"
    echo "-o DIR            Output directory"
    echo
    echo "Other options:"
    echo "-c NUMERIC        Coverage threshold            [default: 0.60]"
    echo "-t DIR            Identity threshold            [default: 0.95]"
    echo "                  NOTE: This is the same as the default on the PlasmidFinder webserver,"
    echo "                        whereas the default for the command-line program is 0.90."
    echo "-d DIR            Directory with PlasmidFinder DB"
    echo "                  default: /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6/share/plasmidfinder-2.1.6/database"
    echo "-a STRING         Other argument(s) to pass to PlasmidFinder"
    echo "-h                Print this help message and exit"
    echo
    echo "Example:          $0 -i my_genome.fa -o results/plasmidfinder"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "PlasmidFinder documentation: https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/"
    echo "PlasmidFinder paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/"
    echo "PlasmidFinder webserver: https://cge.cbs.dtu.dk/services/PlasmidFinder/"
    echo
}

## Option defaults
min_cov=0.60
min_id=0.95

genome_fa=""
db_dir="" && db_arg=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:d:c:t:o:a:h' flag; do
    case "${flag}" in
        i) genome_fa="$OPTARG" ;;
        c) min_cov="$OPTARG" ;;
        t) min_id="$OPTARG" ;;
        d) db_dir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------
## Check input
[[ ! -f "$genome_fa" ]] && echo "## ERROR: Input file (-i) $genome_fa does not exist" >&2 && exit 1

## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6

## Bash strict mode
set -euo pipefail

## PlasmidFinder DB - is at /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6/share/plasmidfinder-2.1.6/database
if [[ "$db_dir" != "" ]]; then
    [[ ! -d "$db_dir" ]] && echo "## ERROR: Database dir (-d) $db_dir does not exist" >&2 && exit 1
    db_arg="-p $db_dir"
fi

## Report
echo "## Starting script plasmidfinder.sh"
date
echo
echo "## Genome FASTA file:                    $genome_fa"
echo "## Output dir:                           $outdir"
[[ "$db_dir" != "" ]] && echo "## Plasmidfinder database dir:           $db_dir"
echo "## Min coverage threshold:               $min_cov"
echo "## Min identity threshold:               $min_id"
echo "## Other arguments for PlasmidFinder:    $more_args"
echo -e "--------------------\n"


# RUN PLASMIDFINDER ------------------------------------------------------------
## Make output dir
mkdir -p "$outdir"

echo "## Now running PlasmidFinder..."
plasmidfinder.py \
    -i "$genome_fa" \
    -o "$outdir" \
    --mincov "$min_cov" \
    --threshold "$min_id" \
    --extented_output \
    $db_arg \
    $more_args

#? Yes, the 'extented_output' option is misspelled by plasmidfinder.py
#? Default PlasmidFinder coverage threshold (--mincov): 0.60
#? Default PlasmidFinder identity threshold (--threshold): 0.90 -- but on the webserver, it's 0.95!


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script plasmidfinder.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo

#? Please run download-db.sh to download the PlasmidFinder database to /fs/project/PAS0471/jelmer/conda/plasmidfinder-2.1.6/share/plasmidfinder-2.1.6/database.
#? If you have a database in custom path, please use plasmidfinder.py with the option -p.