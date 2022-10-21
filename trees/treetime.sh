#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name=treetime
#SBATCH --output=slurm-treetime-%j.out

# FUNCTIONS --------------------------------------------------------------------
## Help function
print_help() {
    echo
    echo "======================================================================"
    echo "$0: Run TreeTime to date a phylogenetic tree"
    echo "======================================================================"
    echo "USAGE:"
    echo "------------------"
    echo "$0 -i <alignment-file> -d <dates-file> -o <output-dir> ..."
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i FILE           Input alignment: FASTA, Phylip, or VCF"
    echo "    -d FILE           Input CSV/TSV with sampling date for each sample"
    echo "                      Should have columns 'node_name' and 'date',"
    echo "                      with dates in %Y-%m-%d format"
    echo "    -o DIR            Output dir"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -t FILE           Tree file in Nexus or Newick format"
    echo "                      If no tree file is provided, TreeTime will infer one"
    echo "    -a STRING         Other argument(s) to pass to TreeTime"
    echo "    -v                Print the version and exit"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "sbatch $0 -i aln.fasta -d data/meta/dates.tsv -o results/treetime"
    echo
    echo "DOCUMENTATION:"
    echo "------------------"
    echo "Documentation:        https://treetime.readthedocs.io/"
    echo "Repo:                 https://github.com/neherlab/treetime"
    echo "Paper:                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5758920/"
    echo
}

## Load software
load_software() {
    module load python/3.6-conda5.2
    source activate /fs/project/PAS0471/jelmer/conda/treetime-0.9.4
    #? Note: The treetime environment also includes FastTree
}

## Print the version
print_version() {
    load_software
    treetime --version
}


# PARSE OPTIONS ----------------------------------------------------------------
## Option defaults
aln=""
dates=""
outdir=""
tree="" && tree_arg=""
clock_rate="" && clock_rate_arg=""
more_args=""

## Parse command-line options
while getopts ':i:d:t:o:c:a:vh' flag; do
    case "${flag}" in
        i) aln="$OPTARG" ;;
        d) dates="$OPTARG" ;;
        t) tree="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        c) clock_rate="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        v) print_version && exit 0 ;;
        h) print_help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ "$aln" = "" ]] && echo "## ERROR: Please specify an alignment with -i" >&2 && exit 1
[[ "$dates" = "" ]] && echo "## ERROR: Please specify a dates file -d" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$aln" ]] && echo "## ERROR: Input file $aln does not exist" >&2 && exit 1
[[ ! -f "$dates" ]] && echo "## ERROR: Input file $dates does not exist" >&2 && exit 1
[[ "$tree" != "" && ! -f "$tree" ]] && echo "## ERROR: Input file $tree does not exist" >&2 && exit 1

## Build tree argument
[[ $tree != "" ]] && tree_arg="--tree $tree"
[[ $clock_rate != "" ]] && clock_rate_arg="--clock-rate $clock_rate"

## Report
echo "=========================================================================="
echo "## Starting script treetime.sh"
date
echo "=========================================================================="
echo
echo "## Alignment (or VCF) input file:             $aln"
echo "## CSV/TSV dates input file:                  $dates"
echo "## Output dir:                                $outdir"
[[ "$tree" != "" ]] && echo "## Tree input file:                           $tree"
[[ $more_args != "" ]] && echo "## Other arguments for treetime:              $more_args"
echo
echo "## Listing the input files:"
ls -lh "$aln" "$dates"
[[ "$tree" != "" ]] && ls -lh "$tree"
echo "=========================================================================="



# MAIN -------------------------------------------------------------------------
## Create output dir
mkdir -p "$outdir"/logs

## Run treetime
echo "## Now running treetime..."
set -o xtrace
treetime \
    --aln "$aln" \
    $tree_arg $clock_rate_arg \
    --dates "$dates" \
    --outdir "$outdir" \
    --confidence
set +o xtrace

echo "## Treetime version used:"
treetime --version | tee "$outdir"/logs/version.txt


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo
tree "$outdir"
echo
echo "## Show contents of dates.tsv:"
cat -n "$outdir"/dates.tsv
echo 
echo "## Done with script treetime.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
