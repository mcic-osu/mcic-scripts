#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
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
    echo "    -i FILE        Gzipped (R1) FASTQ input file (if paired-end, R2 file name will be inferred)"
    echo "    -o DIR         Output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -q INTEGER     Quality trimming threshold         [default: 20 (also the TrimGalore default)]"
    echo "    -l INTEGER     Minimum read length                [default: 20 (also the TrimGalore default)]"
    echo "    -s             Input is single-end                [default: paired-end]"
    echo "    -f             Don't run FastQC after trimming    [default: run FastQC after trimming]"
    echo
    echo "UTTILITY OPTIONS:"
    echo "    -h             Print this help message and exit"
    echo "    -N             Dry run: don't execute commands, only parse arguments and report"
    echo "    -x             Run the script in debug mode (print all code)"
    echo "    -v             Print the version of TrimGalore and exit"
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
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/trimgalore-0.6.7
}

## Print version
Print_version() {
    module load python/3.6-conda5.2
    trim_galore --version
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
qual=20                 # => 20 is also the TrimGalore default
len=20                  # => 20 is also the TrimGalore default
single_end=false        # => paired-end by default

debug=false
dryrun=false
run_fastqc=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
while getopts ':i:o:q:l:fshNxv' flag; do
    case "${flag}" in
        i) R1_in="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        q) qual="$OPTARG" ;;
        l) len="$OPTARG" ;;
        s) single_end=true ;;
        f) run_fastqc=false ;;
        N) dryrun=true ;;
        x) debug=true ;;
        v) Print_version; exit 0 ;;
        h) Print_help; exit 0 ;;
        \?) Die "Invalid option -$OPTARG" ;;
        :) Die "Option -$OPTARG requires an argument" ;;
    esac
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Software
Load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ "$outdir" = "" ]] && Die "Please specify an outdir with -o"
[[ "$R1_in" = "" ]] && Die "Please specify an input R1 FASTQ file with -i"
[[ ! -f "$R1_in" ]] && Die "Input R1 FASTQ file $R1_in does not exist"

## Output dir for TrimGalore logs
outdir_trim="$outdir"/trimmed
outdir_fastqc="$outdir"/fastqc
outdir_logs="$outdir"/logs

## Get number of threads
set +u
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi
set -u

## FastQC arg
if [[ "$run_fastqc" = true ]]; then
    fastqc_arg1="--fastqc --fastqc_args"
    fastqc_arg2="-t $n_threads --outdir $outdir_fastqc"
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
echo "R1 input file:                     $R1_in"
echo "Base output dir:                   $outdir"
echo
[[ "$single_end" != "true" ]] && echo "R2 input file:                     $R2_in"
echo "Sequence quality threshold:        $qual"
echo "Minimum sequence length:           $len"
echo "Sequences are single-end:          $single_end"
echo "Run FastQC:                        $run_fastqc"
echo
echo "Sample ID:                         $sample_id"
echo "Output dir - FastQC:               $outdir_fastqc"
echo "R1 output file:                    $R1_out"
[[ "$single_end" != "true" ]] && echo "R2 output file:                    $R2_out"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    ## Make output dirs
    mkdir -p "$outdir_trim" "$outdir_fastqc" "$outdir_logs"

    ## Run Trim-Galore
    trim_galore \
        --output_dir "$outdir_trim" \
        --quality "$qual" \
        --length "$len" \
        --gzip \
        -j "$n_threads" \
        $fastqc_arg1 "$fastqc_arg2" $input_arg

    ## Move output files
    echo -e "\n## Listing original output files:"
    ls -lh "$outdir_trim"/"$sample_id"*
    
    echo -e "\n## Moving output files..."
    mv -v "$outdir_trim"/"$(basename "$R1_in")"_trimming_report.txt "$outdir_logs"

    if [ "$single_end" != "true" ]; then
        mv -v "$outdir_trim"/"$sample_id"*_val_1.fq.gz "$R1_out"
        mv -v "$outdir_trim"/"$sample_id"*_val_2.fq.gz "$R2_out"
    else
        R1_basename="$(basename "$R1_in" "$extension")"
        mv -v "$outdir_trim"/"$R1_basename"_trimmed.fq.gz "$R1_out"
    fi
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing FASTQ output files:"
    ls -lh "$R1_out"
    [[ "$single_end" != "true" ]] && ls -lh "$R2_out"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
