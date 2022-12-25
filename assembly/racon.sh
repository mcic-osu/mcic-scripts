#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=170G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=racon
#SBATCH --output=slurm-racon-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         RUN MINIMAP => RACON TO POLISH A GENOME ASSEMBLY"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --assembly <assembly-file> --reads <FASTQ-file> -o <output-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --assembly          <file>  Input assembly: FASTA file (to be corrected)"
    echo "  --reads             <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --iterations        <int>   Number of Racon iterations (1 or 2)     [default: 2]"
    echo "  --minimap_preset    <str>   Minimap preset                          [default: 'map-ont']"
    echo "  --more_args_racon   <str>   Quoted string with additional argument(s) to pass to Racon"
    echo "  --more_args_minimap <str>   Quoted string with additional argument(s) to pass to Minimap"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v/--version                Print the version of Racon and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/my.fastq -a results/my.bam -r results/assembly.fasta -o results/racon"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Racon: https://github.com/lbcb-sci/racon"
    echo
}

# Load software
Load_racon() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/racon-1.5.0
}

Load_minimap() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/minimap2-2.24
}

# Print version
Print_version() {
    echo "# Racon:"
    Load_racon
    racon --version
    echo "# Minimap:"
    Load_minimap
    minimap2 --version
}

# Print help for the focal program
Print_help_program() {
    Load_software
    racon --help
}

# Print SLURM job resource usage info
Resource_usage() {
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
}

# Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Job ID:                       $SLURM_JOB_ID"
    echo "Job name:                     $SLURM_JOB_NAME"
    echo "Memory in MB (per node):      $SLURM_MEM_PER_NODE"
    echo "CPUs per task:                $SLURM_CPUS_PER_TASK"
    echo "Nr of tasks:                  $SLURM_NTASKS"
    echo "Account (project):            $SLURM_JOB_ACCOUNT"
    echo "Time limit:                   $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

# Set the number of threads/CPUs
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

# Recource usage information
Time() {
/usr/bin/time -f \
    '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Avg Mem: %t K    Exit status: %x \n' \
    "$@"
}   

# Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}

# Run Racon
Run_racon() {
    assembly_in=${1:-none}
    align=${2:-none}
    assembly_out=${3:-none}

    [[ $assembly_in = "none" ]] && Die "No assembly for function Run_racon"
    [[ $align = "none" ]] && Die "No alignments for function Run_racon"
    [[ $assembly_out = "none" ]] && Die "No outfile for function Run_racon"

    Load_racon

    ${e}Time racon \
        "$reads" \
        "$align" \
        "$assembly_in" \
        -t "$threads" \
        $more_args_racon > "$assembly_out"

    echo
    date
}

# Run Minimap
Run_minimap() {
    assembly=${1:-none}
    align_out=${2:-none}

    [[ $assembly = "none" ]] && Die "No assembly for function Run_minimap"
    [[ $align_out = "none" ]] && Die "No outfile for function Run_minimap"

    Load_minimap
    
    ${e}Time minimap2 \
        -x "$minimap_preset" \
        -t "$threads" \
        -a \
        $more_args_minimap \
        "$assembly" \
        "$reads" \
        > "$align_out"
    
    #? NOTE: Using '--threads' instead of '-t' doesn't seem to work! => Only uses 1 thread

    echo
    date
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
minimap_preset="map-ont"
iterations=2

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
reads=""
assembly_in=""
more_args_racon=""
more_args_minimap=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --reads )               shift && reads=$1 ;;
        --assembly )            shift && assembly_in=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --minimap_preset )      shift && minimap_preset=$1 ;;
        --iterations )          shift && iterations=$1 ;;
        --more_args_racon )     shift && more_args_racon=$1 ;;
        --more_args_minimap )   shift && more_args_minimap=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true && e="echo ";;
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Set nr of threads
Set_threads

# Bash script settings
set -euo pipefail

# Define output file
assembly_ext=$(echo "$assembly_in" | sed -E 's/.*(\.fn?a?s?t?a$)/\1/')
assembly_id=$(basename "$assembly_in" "$assembly_ext")
assembly_out1="$outdir"/"$assembly_id"_racon1.fasta
[[ "$iterations" = "2" ]] && assembly_out2="$outdir"/"$assembly_id"_racon2.fasta
align_1="$outdir"/minimap/"$assembly_id"_iter1.sam
[[ "$iterations" = "2" ]] && align_2="$outdir"/minimap/"$assembly_id"_iter2.sam

# Check input
[[ $reads = "" ]] && Print_args "$all_args" && Die "Please specify a file with input reads with -i"
[[ $assembly_in = "" ]] && Print_args "$all_args" && Die "Please specify an input assembly with -r"
[[ $outdir = "" ]] && Print_args "$all_args" && Die "Please specify an output dir with -o"
[[ ! -f $reads ]] && Die "Input FASTQ file $reads does not exist"
[[ ! -f $assembly_in ]] && Die "Input assembly file $assembly_in does not exist"
[[ "$iterations" != 1 && "$iterations" != 2 ]] && Die "Number of iterations should be 1 or 2 (You asked for $iterations)" 

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RACON.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input reads (FASTQ) file:         $reads"
echo "Input assembly (FASTA) file:      $assembly_in"
echo "Output dir:                       $outdir"
echo "Nr of Racon iterations:           $iterations"
echo "Minimap preset:                   $minimap_preset"
[[ $more_args_racon != "" ]] && echo "Other arguments for Racon:    $more_args_racon"
[[ $more_args_minimap != "" ]] && echo "Other arguments for Minimap:  $more_args_minimap"
echo "Assembly after Racon iteration 1: $assembly_out1"
[[ "$iterations" = 2 ]] && echo "Assembly after Racon iteration 2: $assembly_out2"
echo "SAM after Minimap iteration 1:    $align_1"
[[ "$iterations" = 2 ]] && echo "SAM after Minimap iteration 2:    $align_2"
echo
echo "Listing the input files:"
ls -lh "$reads" "$assembly_in" 
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
${e}mkdir -pv "$outdir"/logs "$outdir"/minimap

# Minimap iteration 1
if [[ ! -s "$align_1" ]]; then
    echo -e "\n# Now running the first iteration of Minimap..."
    ${e}Run_minimap "$assembly_in" "$align_1"
else
    echo -e "\n# Minimap SAM from iteration 1 exists, skipping step..."
    ls -lh "$align_1"
fi

# Racon iteration 1
if [[ ! -s "$assembly_out1" ]]; then
    echo -e "\n====================================================================="
    echo -e "# Now running the first iteration of Racon..."
    ${e}Run_racon "$assembly_in" "$align_1" "$assembly_out1"
else
    echo -e "\n# Assembly from Racon iteration 1 exists, skipping step..."
    ls -lh "$assembly_out1"
fi

if [[ "$iterations" -eq 2 ]]; then
    # Minimap iteration 2
    if [[ ! -s "$align_2" ]]; then
        echo -e "\n====================================================================="
        echo -e "\n# Now running the second iteration of Minimap..."
        ${e}Run_minimap "$assembly_out1" "$align_2"
    else
        echo -e "\n# Minimap SAM from iteration 2 exists, skipping step..."
        ls -lh "$align_2"
    fi

    # Racon iteration 2
    if [[ ! -s "$assembly_out2" ]]; then
    echo -e "\n====================================================================="
    echo "# Now running the second iteration of Racon..."
        ${e}Run_racon "$assembly_out1" "$align_2" "$assembly_out2"
    else
        echo -e "\n# Assembly from Racon iteration 2 exists, skipping step..."
        ls -lh "$assembly_out2"
    fi
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    echo
    [[ "$slurm" = true ]] && Resource_usage
    echo
fi
echo "# Done with script"
date
