#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=evigene
#SBATCH --output=slurm-evigene-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  TODO FUNCTION OF THIS SCRIPT"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  FASTA file containing multiple, concatenated transcriptome assemblies"
    echo "                          NOTE: Concatenate multiple assemblies FASTAs into a single one prior to running this script."
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "                          NOTE: Any files present in this dir will be removed prior to running EviGene!"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Evigene"
    echo "  --min_cds_length <int>  Minimum CDS length in bp                    [default: 350]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Evigene and exit"
    echo
    echo "OUTPUT:"
    echo "  - The final full transcriptome file will be '<outdir>final/evigene_all.fasta'"
    echo "  - The equivalent file with only primary transcript will be '<outdir>/final/evigene_primarytrans.fasta'"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i result/concat_assembly.fasta -o results/evigene"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - EviGene documentation: http://arthropods.eugenes.org/EvidentialGene/evigene/"
    echo "  - EviGene output:        http://arthropods.eugenes.org/EvidentialGene/evigene/docs/EvigeneR/evigene4_outputs_brief.txt"
    echo "  - EviGene update note:   http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_update2020march.html"
    echo "  - EviGene paper:         https://www.biorxiv.org/content/10.1101/829184v1"
    echo "  - Evigene repo:          https://sourceforge.net/projects/evidentialgene"
    echo
}


## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    CONDA_ENV_DIR=/fs/project/PAS0471/jelmer/conda/evigene
    source activate "$CONDA_ENV_DIR"
    EVIGENE="$CONDA_ENV_DIR"/bin/prot/tr2aacds.pl
}

## Print help for the focal program
Print_help_program() {
    Load_software
    "$EVIGENE"
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

## Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Account (project):    $SLURM_JOB_ACCOUNT"
    echo "Job ID:               $SLURM_JOB_ID"
    echo "Job name:             $SLURM_JOB_NAME"
    echo "Memory (per node):    $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

## Set the number of threads/CPUs
Set_threads() {
    set +u
    if [[ "$slurm" = true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            threads="$SLURM_NTASKS"
        else 
            echo "WARNING: Can't detect nr of threads, setting to 1"
            threads=1
        fi
    else
        threads=1
    fi
    set -u
}

## Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo >&2
    echo "=====================================================================" >&2
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option" >&2
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'" >&2
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    echo -e "\nEXITING..." >&2
    echo "=====================================================================" >&2
    echo >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
min_cds=350

mem=4000               # In MB; only applies if not a Slurm job
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
outdir=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --min_cds_length )  shift && min_cds_length=$1 ;;
        --more-args )       shift && more_args=$1 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads
[[ "$slurm" = true ]] && mem="$SLURM_MEM_PER_NODE"

## Bash script settings
set -euo pipefail

## Define output files (note: this is within the $outdir, where we will cd to)
assembly_all=final/evigene_all.fasta
assembly_primary=final/evigene_primarytrans.fasta

## Check input
[[ "$infile" = "" ]] && Die "Please specify an input file with -i/--infile" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && Die "Input file $infile does not exist"

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT EVIGENE.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $infile"
echo "Output dir:                       $outdir"
echo "Minimum CDS size:                 $min_cds_length"
[[ $more_args != "" ]] && echo "Other arguments for Evigene:      $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Remove output dir if it already exists (EviGene behaves weirdly when there are files present!)
    [[ "$outdir" != "" && -d "$outdir" ]] && rm -r "${outdir:?}"

    ## Create the output directory
    ${e}mkdir -p "$outdir"/logs "$outdir"/final

    ## Copy input file to outdir
    echo "## Copying input FASTA to output dir..."
    infile_base=$(basename "$infile")
    cp -v "$infile" "$outdir"

    ## Move into the output dir
    cd "$outdir" || exit

    echo -e "\n## Now running Evigene..."
    ${e}Time $EVIGENE \
        -mrnaseq "$infile_base" \
        -debug \
        -MINCDS "$min_cds" \
        -NCPU "$threads" \
        -MAXMEM "$mem" \
        -log \
        -tidyup \
        $more_args

    ## Copy the final assembly file
    cp -v okayset/merged.okay.mrna "$assembly_all"

    ## Make a separate file with primary transcripts
    awk -v RS='>' '/t1 type/ { print ">" $0 }' "$assembly_all" > "$assembly_primary"

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
