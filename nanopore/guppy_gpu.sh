#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --gpus-per-node=2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=guppy
#SBATCH --output=slurm-guppy_gpu-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "      PERFORM FAST5 => FASTQ BASECALLING WITH GUPPY USING GPUS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> --config <config file> [...]"
    echo "  bash $0 --help"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input FAST5 file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  --config        <str>   Config file name"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --barcode_kit   <str>   Barcode set                                 [default: none - no demultiplexing]"
    echo "  --min_qscore    <int>   Minimum quality-score                       [default: whatever is specified in the focal config file]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Guppy"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Guppy and exit"
    echo "  -v/--version            Print the version of Guppy and exit"
    echo
    echo "HARDCODED GUPPY PARAMETERS:"
    echo "  - --chunks_per_runner=256"
    echo "  - --records_per_fastq=0"
    echo "  - --compress_fastq"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/minion/my.fast5 -o results/guppy --config dna_r9.4.1_450bps_sup.cfg"
    echo
    echo "This script will run Guppy for one FAST5 file, so you should loop over the FAST5 files, e.g.:"
    echo "  for fast5 in data/minion/fast5/*fast5; do"
    echo '      outdir=results/guppy/$(basename "$infile" .fast5)'
    echo '      sbatch $0 -i $fast5 -o results/guppy --config dna_r9.4.1_450bps_sup.cfg'
    echo "  done"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revan_14dec2018/run-guppy-on-linux"
    echo "  - https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revan_14dec2018/setting-up-a-run-configurations-and-parameters"
    echo
}

# Print version
Print_version() {
    set +e
    $GUPPY_BIN --version | head -1
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    $GUPPY_BIN --help
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

# Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Account (project):    $SLURM_JOB_ACCOUNT"
    echo "Job ID:               $SLURM_JOB_ID"
    echo "Job name:             $SLURM_JOB_NAME"
    echo "GPUs per node:        $SLURM_GPUS_PER_NODE"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

# Set the number of threads/CPUs
Set_threads() {
    set +u
    if [[ "$slurm" = true ]]; then
        if [[ -n "$SLURM_GPUS_PER_NODE" ]]; then
            threads="$SLURM_GPUS_PER_NODE"
        else 
            echo "WARNING: Can't detect nr of GPUs, setting to 1"
            threads=1
        fi
    else
        threads=1
    fi
    set -u
}

# Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
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
# Software
GUPPY_DIR=/fs/project/PAS0471/jelmer/software/guppy-6.4.2 # https://community.nanoporetech.com/downloads
GUPPY_BIN="$GUPPY_DIR"/bin/guppy_basecaller
module load cuda

# Hardcoded parameters
chunks_per_runner=256
records_per_fastq=0

# Option defaults
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
infile=""
outdir=""
config=""
min_qscore="" && qscore_arg=""
barcode_kit="" && barcode_arg=""
more_args=""

# Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --config )          shift && config=$1 ;;
        --barcode_kit )     shift && barcode_kit=$1 ;;
        --min_qscore )      shift && min_qscore=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
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
# Bash script settings
set -euo pipefail

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Set nr of threads
Set_threads

# Build input argument
indir=$(dirname "$infile")
infile_base=$(basename "$infile")
fofn="$outdir"/tmp/"$infile_base".fofn
infile_arg="--input_file_list $fofn"

# Optional parameters
[[ $barcode_kit != "" ]] && barcode_arg="--barcode_kits $barcode_kit --trim_barcodes"
[[ $min_qscore != "" ]] && qscore_arg="--qscore_arg $min_qscore"

# Check input
[[ "$infile" = "" ]] && Die "Please specify an input FAST5 file with -i/--infile" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ "$config" = "" ]] && Die "Please specify an input config name with --config" "$all_args"
[[ ! -f "$infile" ]] && Die "Input FAST5 file $infile does not exist"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT GUPPY_GPU.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $infile"
echo "Output dir:                       $outdir"
echo "ONT config name:                  $config"
[[ $min_qscore != "" ]] && echo "Minimum qual score to PASS:       $min_qscore"
[[ $barcode_kit != "" ]] && echo "Barcode kit name:                 $barcode_kit"
[[ $more_args != "" ]] && echo "Other arguments for Guppy:        $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$infile"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory if it doesn't already exist
${e}mkdir -p "$outdir"/tmp "$outdir"/logs

# Create fofn
echo "$infile_base" > "$fofn"

# Run Guppy
echo "Now running Guppy..."
${e}Time $GUPPY_BIN \
    --input_path "$indir" \
    $infile_arg \
    --save_path "$outdir" \
    --config "$config" \
    --compress_fastq \
    --records_per_fastq "$records_per_fastq" \
    --gpu_runners_per_device "$threads" \
    --as_gpu_runners_per_device "$threads" \
    --chunks_per_runner "$chunks_per_runner" \
    --num_callers "$threads" \
    --device "cuda:0" \
    $qscore_arg \
    $barcode_arg \
    $more_args

#? GPU parameters
# - https://github.com/colindaven/guppy_on_slurm/blob/master/runbatch_gpu_guppy.sh
#   gpu_params='--device "cuda:0" --num_callers 4 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 3000'
# - See "Tuning GPU parameters for Guppy performance" https://hackmd.io/@Miles/S12SKP115

#? To print config names for flowcell + kit combs:
#$ guppy_basecaller --print_workflows
# - For kit LSK109 and flowcell FLO-MIN106, this is dna_r9.4.1_450bps_hac
# - See also the available config files is /fs/ess/PAS0471/jelmer/software/guppy-6.4.2/data/


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
    echo -e "\n# Listing the 'pass' FASTQ file:"
    ls -lh "$PWD"/"$outdir"/pass
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
