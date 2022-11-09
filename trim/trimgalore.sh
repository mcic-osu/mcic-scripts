#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm-trimgalore-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help
Print_help() {
    echo
    echo "======================================================================"
    echo "                             $0"
    echo "     Run TrimGalore for one (single-end) or a pair of FASTQ files"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-FASTQ> -o <outdir> [...]"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Gzipped (R1) FASTQ input file (if paired-end, R2 file name will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -q/--quality    <int>   Quality trimming threshold         [default: 20 (also the TrimGalore default)]"
    echo "  -l/--length     <int>   Minimum read length                [default: 20 (also the TrimGalore default)]"
    echo "  -s/--single_end         Input is single-end                [default: paired-end]"
    echo "  --no_fastqc             Don't run FastQC after trimming    [default: run FastQC after trimming]"
    echo "  --more_args             Additional arguments to pass to TrimGalore"
    echo
    echo "UTTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for TrimGalore and exit"
    echo "  -v/--version            Print the version of TrimGalore and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq/S01_R1.fastq.gz -o results/trimgalore"
    echo
    echo "OUTPUT:"
    echo "  Within the specified output dir (-o), results will be in 3 subdirectories:"
    echo "    - Trimmed FASTQ files will be in '<outdir>/trimmed/'"
    echo "    - FastQC output be in '<outdir>/fastqc/'"
    echo "    - Trimmomatic log files will be in '<outdir>/logs/'"
    echo
    echo "DOCUMENTATION:"
    echo "  - https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md"
    echo
}

## Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/trimgalore-0.6.7
    set -u
}

## Print version
Print_version() {
    Load_software
    trim_galore --version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    trim_galore --help
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

## Print SLURM job resource usage info
Resource_usage() {
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
}

## Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Job ID:                       $SLURM_JOB_ID"
    echo "Job name:                     $SLURM_JOB_NAME"
    echo "Memory in GB (per node):      $SLURM_MEM_PER_NODE"
    echo "CPUs per task:                $SLURM_CPUS_PER_TASK"
    echo "Nr of tasks:                  $SLURM_NTASKS"
    echo "Account (project):            $SLURM_JOB_ACCOUNT"
    echo "Time limit:                   $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' / '--help' option"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}

## Recource usage information
Time() {
/usr/bin/time -f \
    '\n# Ran the command:\n%C \n# Run stats:\nTime: %E   CPU: %P    Max mem: %M K    Avg Mem: %t K    Exit status: %x \n' \
    "$@"
}   


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
quality=20                 # => 20 is also the TrimGalore default
length=20                  # => 20 is also the TrimGalore default
single_end=false        # => paired-end by default

debug=false
dryrun=false && e=""
run_fastqc=true
slurm=true

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
R1_in=""
outdir=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )             shift && R1_in=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -q | --quality )        shift && quality=$1 ;;
        -l | --length )         shift && length=$1 ;;
        -s | --single_end )     single_end=true ;;
        --no_fastqc )           run_fastqc=false ;;
        --more_args )           shift && more_args=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true && e="echo ";;
        -v | --version )        Print_version; exit 0;;
        -h )                    Print_help; exit 0;;
        --help )                Print_help_program; exit 0;;
        * )                     Die "Invalid option $1" "$all_args";;
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

## Bash strict settings
set -euo pipefail

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

## Check input
[[ "$outdir" = "" ]] && Die "Please specify an outdir with -o"
[[ "$R1_in" = "" ]] && Die "Please specify an input R1 FASTQ file with -i"
[[ ! -f "$R1_in" ]] && Die "Input R1 FASTQ file $R1_in does not exist"

## Output dir for TrimGalore logs
outdir_trim="$outdir"/trimmed
outdir_fastqc="$outdir"/fastqc
outdir_logs="$outdir"/logs

## FastQC arg
if [[ "$run_fastqc" = true ]]; then
    fastqc_arg1="--fastqc --fastqc_args"
    fastqc_arg2="-t $threads --outdir $outdir_fastqc"
else
    fastqc_arg1=""
    fastqc_arg2=""
fi

## Get file extension (.fastq.gz or .fq.gz)
extension=$(echo "$R1_in" | sed -E 's/.*(\.fa?s?t?q\.gz$)/\1/')

## Get R2 file, create input argument, define output files
if [ "$single_end" != "true" ]; then

    # Paired-end sequences
    R1_suffix=$(echo "$R1_in" | sed -E "s/.*(_R?1)_?[[:digit:]]*$extension/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    input_arg="--paired $R1_in $R2_in"
    
    [[ ! -f "$R2_in" ]] && Die "Input R2 FASTQ file $R2_in does not exist"
    [[ "$R1_in" = "$R2_in" ]] && Die "Input R1 and R2 FASTQ files are the same file"

    sample_id=$(basename "$R1_in" | sed -E "s/${R1_suffix}_?[[:digit:]]*${extension}//")
    R1_out="$outdir_trim"/"$sample_id"_R1.fastq.gz
    R2_out="$outdir_trim"/"$sample_id"_R2.fastq.gz
else
    
    # Single-end sequences
    input_arg="$R1_in"
    R1_suffix=""

    sample_id=$(basename "$R1_in" | sed "s/${R1_suffix}${extension}//")
    R1_out="$outdir_trim"/"$sample_id".fastq.gz
fi

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT TRIMGALORE.SH"
date
echo "=========================================================================="
echo "All arguments:                    $all_args"
echo "R1 input file:                    $R1_in"
echo "Base output dir:                  $outdir"
echo "Sequence quality threshold:       $quality"
echo "Minimum sequence length:          $length"
echo "Sequences are single-end:         $single_end"
echo "Run FastQC:                       $run_fastqc"
echo
[[ $more_args != "" ]] && echo "Other arguments for TrimGalore:   $more_args"
[[ "$single_end" != "true" ]] && echo "R2 input file (inferred):         $R2_in"
echo "Sample ID (inferred):             $sample_id"
echo "Output dir - FastQC:              $outdir_fastqc"
echo "R1 output file:                   $R1_out"
[[ "$single_end" != "true" ]] && echo "R2 output file:                    $R2_out"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Make output dirs
${e}mkdir -p "$outdir_trim" "$outdir_fastqc" "$outdir_logs"

## Run Trim-Galore
${e}Time \
    trim_galore \
        --output_dir "$outdir_trim" \
        --quality "$quality" \
        --length "$length" \
        --gzip \
        -j "$threads" \
        $more_args \
        $fastqc_arg1 "$fastqc_arg2" $input_arg
    
## Move output files
echo -e "\n# Moving output files..."
${e}mv -v "$outdir_trim"/"$(basename "$R1_in")"_trimming_report.txt "$outdir_logs"

if [ "$single_end" != "true" ]; then
    ${e}mv -v "$outdir_trim"/"$sample_id"*_val_1.fq.gz "$R1_out"
    ${e}mv -v "$outdir_trim"/"$sample_id"*_val_2.fq.gz "$R2_out"
else
    R1_basename="$(basename "$R1_in" "$extension")"
    ${e}mv -v "$outdir_trim"/"$R1_basename"_trimmed.fq.gz "$R1_out"
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing FASTQ output files:"
    ls -lh "$R1_out"
    [[ "$single_end" != "true" ]] && ls -lh "$R2_out"
    echo
    [[ "$slurm" = true ]] && Resource_usage
fi
echo
echo "# Done with script"
date
