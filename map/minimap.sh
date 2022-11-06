#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --job-name=minimap
#SBATCH --output=slurm-minimap-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "             MAP READS TO A GENOME WITH MINIMAP2"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-dir> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--reads      <file>      Input FASTQ file"
    echo "  -r/--reference  <file>      Input reference"
    echo "  -o/--outdir     <dir>       Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -O/--outfile_type <string>  Output file type: 'sam' or 'paf'        [default: 'sam']"
    echo "  -x/--preset     <string>    Preset: read and operation type, see below for list [default: 'map-ont']" 
    echo "  -a/--more_args  <string>    Quoted string with additional argument(s) to pass to Minimap2"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -v/--version                Print the version of Minimap2 and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i TODO -o results/TODO "
    echo "  sbatch $0 -i TODO -o results/TODO -a \"-x TODO\""
    echo
    echo "MINIMAP2 PRESET OPTIONS:"
    echo "  - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping"
    echo "  - map-hifi - PacBio HiFi reads vs reference mapping"
    echo "  - ava-pb/ava-ont - PacBio/Nanopore read overlap"
    echo "  - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence"
    echo "  - splice/splice:hq - long-read/Pacbio-CCS spliced alignment"
    echo "  - sr - genomic short-read mapping"
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
    source activate /fs/ess/PAS0471/jelmer/conda/minimap2-2.24
}

## Print version
Print_version() {
    Load_software
    minimap2 --version
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
outfile_type=sam && outfile_arg="-a"
preset="map-ont"

debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
reads=""
reference=""
outdir=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )          shift && reads=$1 ;;
        -r | --reference )      shift && reference=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -O | --outfile_type )   shift && outfile_type=$1 ;;
        -x | --preset )         shift && preset=$1 ;;
        -a | --more_args )      shift && more_args=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true ;;
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

## Output file option
[[ "$outfile_type" = paf ]] && outfile_arg=""

## Define output file
reference_ext=$(echo "$reference" | sed -E 's/.*(\.fn?a?s?t?a$)/\1/')
outfile="$outdir"/$(basename "$reference" "$reference_ext")."$outfile_type"

## Check input
[[ $reference = "" ]] && Die "Please specify an input reference genome FASTA with -i"
[[ $reads = "" ]] && Die "Please specify an input FASTQ file with -I"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $reference ]] && Die "Input reference genome FASTA file $reference does not exist"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT MINIMAP.SH"
date
echo "=========================================================================="
echo "Input reads:                 $reads"
echo "Input reference:             $reference"
echo "Output dir:                  $outdir"
echo "Output file type:            $outfile_type"
echo "Minimap2 preset:             $preset"
echo
echo "Output file:                 $outfile"
[[ $more_args != "" ]] && echo "Other arguments for Minimap2:    $more_args"
echo "Listing input files:"
ls -lh "$reads" "$reference"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN\n"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    
    ## Create the output directory
    mkdir -p "$outdir"/logs

    ## Run
    echo -e "\n## Now running Minimap2..."
    [[ "$debug" = false ]] && set -o xtrace
    
    minimap2 \
        -x "$preset" \
        -t "$threads" \
        "$outfile_arg" \
        $more_args \
        "$reference" \
        "$reads" \
        > "$outfile"

    # | samtools sort -o "$bam_out" -T reads.tmp -

    [[ "$debug" = false ]] && set +o xtrace

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    echo
    echo "========================================================================="
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
