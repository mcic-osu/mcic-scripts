#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
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
    echo "         RUN MINIMAP => RACON TO POLISH A GENOME ASSEMBLY"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <FASTQ-file> -r <assembly-file> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--reads      <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  -r/--assembly   <file>  Input assembly: FASTA file (to be corrected)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -x/--minimap_preset <str>  Minimap preset                          [default: 'map-ont']"
    echo "  --more_args_racon  <str>   Quoted string with additional argument(s) to pass to Racon"
    echo "  --more_args_minimap  <str> Quoted string with additional argument(s) to pass to Minimap"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of Racon and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/my.fastq -a results/my.bam -r results/assembly.fasta -o results/racon"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Racon: https://github.com/lbcb-sci/racon"
    echo
}

## Load software
Load_racon() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/racon-1.5.0
}

Load_minimap() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/minimap2-2.24
}

## Print args
Print_args() {
    echo -e "\n# Arguments passed to the script:"
    echo "$*"
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

Run_racon() {
    assembly_in=${1:-none}
    align=${2:-none}
    assembly_out=${3:-none}

    [[ $assembly_in = "none" ]] && Die "No assembly for function Run_racon"
    [[ $align = "none" ]] && Die "No alignments for function Run_racon"
    [[ $assembly_out = "none" ]] && Die "No outfile for function Run_racon"

    Load_racon

    racon \
        "$reads" \
        "$align" \
        "$assembly_in" \
        --threads "$threads" \
        $more_args_racon \
        > "$assembly_out"
}

Run_minimap() {
    assembly=${1:-none}
    align_out=${2:-none}

    [[ $assembly = "none" ]] && Die "No assembly for function Run_minimap"
    [[ $align_out = "none" ]] && Die "No outfile for function Run_minimap"

    Load_minimap
    
    minimap2 \
        -x "$minimap_preset" \
        -t "$threads" \
        -a \
        $more_args_minimap \
        "$assembly" \
        "$reads" \
        > "$align_out"
}

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
minimap_preset="map-ont"
iterations=2

debug=false
dryrun=false && e=""


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
reads=""
assembly_in=""
more_args_racon=""
more_args_minimap=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )          shift && reads=$1 ;;
        -r | --assembly )       shift && assembly_in=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --minimap_preset )      shift && minimap_preset=$1 ;;
        --iterations )          shift && iterations=$1 ;;
        --more_args_racon )     shift && more_args_racon=$1 ;;
        --more_args_minimap )   shift && more_args_minimap=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true && e="echo";;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_args "$all_args"; Die "Invalid option $1" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

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
assembly_ext=$(echo "$assembly_in" | sed -E 's/.*(\.fn?a?s?t?a$)/\1/')
assembly_id=$(basename "$assembly_in" "$assembly_ext")
assembly_out1="$outdir"/"$assembly_id"_racon1.fasta
[[ "$iterations" = "2" ]] && assembly_out2="$outdir"/"$assembly_id"_racon2.fasta
align_1="$outdir"/minimap/"$assembly_id"_iter1.sam
[[ "$iterations" = "2" ]] && align_2="$outdir"/minimap/"$assembly_id"_iter2.sam

## Check input
[[ $reads = "" ]] && Print_args "$all_args" && Die "Please specify a file with input reads with -i"
[[ $assembly_in = "" ]] && Print_args "$all_args" && Die "Please specify an input assembly with -r"
[[ $outdir = "" ]] && Print_args "$all_args" && Die "Please specify an output dir with -o"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"
[[ ! -f $assembly_in ]] && Die "Input assembly file $assembly_in does not exist"
[[ "$iterations" != 1 && "$iterations" != 2 ]] && Die "Number of iterations should be 1 or 2 (You asked for $iterations)" 

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RACON.SH"
date
echo "=========================================================================="
echo "Input reads (FASTQ) file:         $reads"
echo "Input assembly (FASTA) file:      $assembly_in"
echo "Output dir:                       $outdir"
echo "Nr of Racon iterations:           $iterations"
echo "Minimap preset:                   $minimap_preset"
[[ $more_args_racon != "" ]] && echo "Other arguments for Racon:    $more_args_racon"
[[ $more_args_minimap != "" ]] && echo "Other arguments for Minimap:  $more_args_minimap"
echo
echo "Output file after 1st Racon iteration: $assembly_out1"
[[ "$iterations" = "2" ]] && echo "Output file after 2nd Racon iteration: $assembly_out2"
echo
echo "Listing input files:"
ls -lh "$reads" "$assembly_in" 
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
"${e}"mkdir -p "$outdir"/logs "$outdir"/minimap

## Run
[[ "$dryrun" = false ]] && set -o xtrace

echo -e "\n## Now running the first iteration of Minimap..."
"${e}"Run_minimap "$assembly_in" "$align_1"

echo -e "\n## Now running the first iteration of Racon..."
"${e}"Run_racon "$assembly_in" "$align_1" "$assembly_out1"

if [[ "$iterations" = "2" ]]; then
    echo -e "\n## Now running the second iteration of Minimap..."
    "${e}"Run_minimap "$assembly_out1" "$align_2"

    echo -e "\n## Now running the second iteration of Racon..."
    "${e}"Run_racon "$assembly_out1" "$align_2" "$assembly_out2"
fi

[[ "$debug" = false ]] && set +o xtrace

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "## Version used:"
"${e}"Print_version | tee "$outdir"/logs/version.txt
echo -e "\n## Listing files in the output dir:"
"${e}"ls -lhd "$PWD"/"$outdir"/*
echo
"${e}"sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
echo "## Done with script"
date
