#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=43
#SBATCH --mem=172G
#SBATCH --job-name=redbean
#SBATCH --output=slurm-redbean-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "      RUN REDBEAN (WTDBG2) TO ASSEMBLE A GENOME WITH LONG READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-file> -o <output-dir> -p <genome-ID> -g <genome-size> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i/--infile       FILE        Input file"
    echo "    -o/--outdir       DIR         Output dir (will be created if needed)"
    echo "    -p/--genome_id    STRING      Genome ID (prefix) for output files"
    echo "    -g/--genome_size  STRING      Genome size estimate, e.g '4.6m'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    --kmer-subsampling INTEGER    K-mer subsampling rate, lower nr is less subsampling  [default: 4]"
    echo "                                  Less subsampling may improve the assembly for low cov sequencing, but will take longer"
    echo "    --aln-noskip                  Even if a read was contained in previous alignment, still align it against other reads"
    echo "                                  May improve the assembly for low cov sequencing, but will take longer"
    echo "    --lowcov_edge                 For low cov. sequencing: set min edge depth to 2 and use --rescue-low-cov-edges"
    echo "    -a/--more_args    STRING      Quoted string with additional argument(s) to pass to Redbean"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h/--help                     Print this help message and exit"
    echo "    -N/--dryrun                   Dry run: don't execute commands, only parse arguments and report"
    echo "    -x/--debug                    Run the script in debug mode (print all code)"
    echo "    -v/--version                  Print the version of Redbean and exit"
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
    echo "    - https://github.com/ruanjue/wtdbg2"
    echo "    - https://github.com/ruanjue/wtdbg2/blob/master/README-ori.md"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/wtdbg-2.5
}

## Print version
Print_version() {
    Load_software
    wtdbg2 -V
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
kmer_subsampling=4
aln_noskip=false && aln_skip_arg=""
read_type=ont                #  "rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS reads and "ont" for Oxford Nanopore

debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
outdir=""
genome_size=""
lowcov_edge=false && edge_arg=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -p | --genome_id )      shift && genome_id=$1 ;;
        -g | --genome_size )    shift && genome_size=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        --kmer-subsampling )    shift && kmer_subsampling=$1 ;;
        --aln-noskip )          aln_noskip=true ;;
        --lowcov_edge )         lowcov_edge=true ;;
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

## Build other args
[[ "$aln_noskip" = true ]] && aln_skip_arg="--aln-noskip"
[[ "$lowcov_edge" = true ]] && edge_arg="--edge-min 2 --rescue-low-cov-edges"

## Check input
[[ $infile = "" ]] && Die "Please specify an input file with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT REDBEAN.SH"
date
echo "=========================================================================="
echo "Input file:                  $infile"
echo "Output dir:                  $outdir"
echo "Genome ID (prefix):          $genome_id"
echo "Genome size:                 $genome_size"
[[ $more_args != "" ]] && echo "Other arguments for Redbean:    $more_args"
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

    echo -e "\n## Now running 'wtdbg2'..."
    [[ "$debug" = false ]] && set -o xtrace
    wtdbg2 \
        -x $read_type \
        -g "$genome_size" \
        -t "$n_threads" \
        -i "$infile" \
        -fo "$outdir"/"$genome_id" \
        --kmer-subsampling "$kmer_subsampling" \
        -f \
        $aln_skip_arg $edge_arg \
        $more_args
    
    echo -e "\n## Now running the consensus step..."
    wtpoa-cns \
        -i "$outdir"/"$genome_id".ctg.lay.gz \
        -fo "$outdir"/"$genome_id".ctg.fa \
        -t "$n_threads"

    [[ "$debug" = false ]] && set +o xtrace
fi

## -f forces overwrite

## Low coverage options:
# - By default, wtdbg2 samples a fourth of all k-mers by their hashcodes.
#   For data of relatively low coverage, you may increase this sampling rate by reducing -S.
# - Option -A/--aln-noskip also helps relatively low coverage data at the cost of performance. 
#  - For low coverage: Try --edge-min 2 --rescue-low-cov-edges.


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
