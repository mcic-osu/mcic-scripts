#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=sortmerna
#SBATCH --output=slurm-sortmerna-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo " Run SortMeRNA to sort RNAseq reads into rRNA-derived and other reads"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input R1 FASTQ file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Input R1 FASTQ file (name of R2 will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --repo                  Directory with SortMeRNA repo (for reference FASTA files)  [default: download repo]"
    echo "  --leave-interleaved     Don't 'de-interleave' output FASTQ file                    [default: de-interleave]"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to SortmeRNA"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for SortmeRNA and exit"
    echo "  -v/--version            Print the version of SortmeRNA and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/S1_R1.fastq.gz -o results/sortmerna"
    echo
    echo "OUTPUT:"
    echo "  - Aligned sequences will be placed in <output-dir>/mapped"
    echo "  - Non-aligned sequences will be placed in <output-dir>/unmapped"
    echo "Output sequence files will keep sample identifiers,"
    echo "so the script can be run for multiple samples using the same <output-dir> (-o)"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/biocore/sortmerna"
    echo "  - Paper: https://academic.oup.com/bioinformatics/article/28/24/3211/246053"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/sortmerna-env # NOTE: this env also had BBMap installed
    set -u
}

# Print version
Print_version() {
    set -e
    Load_software
    sortmerna --version
    set +e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    sortmerna --help
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
}

# Print SLURM job requested resources
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
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
deinterleave=true

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
R1=""
outdir=""
repo_dir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --leave-interleaved ) deinterleave=false ;;
        --repo )            repo_dir=false ;;
        --more-args )       shift && more_args=$1 ;;
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

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Infer the name of the R2 file
file_ext=$(basename "$R1" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2=${R1/$R1_suffix/$R2_suffix}
sampleID=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")

# Check input
[[ "$R1" = "" ]] && Die "Please specify an input R1 file with -i/--R1" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$R1" ]] && Die "Input file $R1 does not exist"
[[ ! -f "$R2" ]] && Die "Input file $R2 does not exist"
[[ "$R1" = "$R2" ]] && Die "Input file R1 and R2 refer to the same path ($R1)"

# Define output files
outdir_full="$outdir"/"$sampleID"
out_mapped_raw="$outdir_full"/mapped_raw/"$sampleID"
out_unmapped_raw="$outdir_full"/unmapped_raw/"$sampleID"

R1_mapped="$outdir"/mapped/"$sampleID""$R1_suffix".fastq.gz
R2_mapped="$outdir"/mapped/"$sampleID""$R2_suffix".fastq.gz
R1_unmapped="$outdir"/unmapped/"$sampleID""$R1_suffix".fastq.gz
R2_unmapped="$outdir"/unmapped/"$sampleID""$R2_suffix".fastq.gz

# Reference FASTA files (to be downloaded)
[[ $repo_dir = "" ]] && repo_dir="$outdir"/"$sampleID"/sortmerna_repo
ref_18s="$repo_dir"/data/rRNA_databases/silva-euk-18s-id95.fasta
ref_28s="$repo_dir"/data/rRNA_databases/silva-euk-28s-id98.fasta

# Report
echo
echo "=========================================================================="
echo "                 STARTING SCRIPT SORTMERNA.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "R1 FASTQ file:                    $R1"
echo "Output dir:                       $outdir"
echo
echo "SortMeRNA repo dir:               $repo_dir"
echo "Deinterleave FASTQ files:         $deinterleave"
echo "18S reference file:               $ref_18s"
echo "28S reference file:               $ref_28s"
echo "R2 FASTQ file:                    $R2"
[[ $more_args != "" ]] && echo "Other arguments for SortMeRNA:    $more_args"
echo "Number of threads/cores:          $threads"
echo "Listing the input file(s):"
ls -lh "$R1" "$R2"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    # Make output dirs if needed
    mkdir -p "$outdir"/mapped "$outdir"/unmapped \
        "$outdir_full"/mapped_raw "$outdir_full"/unmapped_raw \
        "$outdir"/logs

    # GET DATABASE FILES -----------------------------------------------------------
    # Clone sortmerna repo to get db FASTA files
    if [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]]; then
        n_seconds=$(( RANDOM % 50 + 1 ))
        sleep "$n_seconds"s # Sleep for a while so git doesn't error when running this multiple times in parallel
        
        mkdir -p "$repo_dir"
        echo "# Cloning sortmerna repo..."
        [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]] && git clone https://github.com/biocore/sortmerna "$repo_dir"
    fi

    # Check that db files are there
    [[ ! -f "$ref_18s" ]] && Die "18s reference FASTA file $ref_18s not found"
    [[ ! -f "$ref_28s" ]] && Die "28s reference FASTA file $ref_28s not found"

    # RUN SortMeRNA ----------------------------------------------------------------
    echo -e "# Starting SortMeRNA run....\n"
    Time sortmerna \
        --ref "$ref_18s" \
        --ref "$ref_28s" \
        --reads "$R1" \
        --reads "$R2" \
        --fastx \
        --aligned "$out_mapped_raw" \
        --other "$out_unmapped_raw" \
        --workdir "$outdir_full" \
        --paired_in \
        --threads "$threads" \
        $more_args

    #?--paired_in Flags the paired-end reads as Aligned, when either of them is Aligned.

    # CONVERTING INTERLEAVED FASTQ BACK TO SEPARATED -------------------------------
    if [[ "$deinterleave" = true ]]; then
        echo -e "\n# Deinterleaving mapped reads..."
        Time reformat.sh \
            in="$out_mapped_raw".fq.gz \
            out1="$R1_mapped" \
            out2="$R2_mapped"

        echo -e "\n# Deinterleaving unmapped reads..."
        Time reformat.sh \
            in="$out_unmapped_raw".fq.gz \
            out1="$R1_unmapped" \
            out2="$R2_unmapped"
        
        echo
    else
        mv -v "$out_mapped_raw".fq.gz "$outdir"/mapped
        mv -v "$out_unmapped_raw".fq.gz "$outdir"/unmapped
    fi

    # HOUSEKEEPING -----------------------------------------------------------------
    # Move log files to main dir
    mv "$outdir_full"/mapped_raw/"$sampleID"*log "$outdir"/logs/

    # Remove temporary files
    rm -rv "$outdir_full"/mapped_raw "$outdir_full"/unmapped_raw

    # QUANTIFY MAPPING SUCCESS -----------------------------------------------------
    n_mapped=$(zcat "$R1_mapped" | awk '{ s++ } END{ print s/4 }')
    n_unmapped=$(zcat "$R1_unmapped" | awk '{ s++ } END{ print s/4 }')
    pct=$(python3 -c "print(round($n_mapped / ($n_unmapped + $n_mapped) * 100, 2))")
    echo -e "\nNumber of reads mapped/unmapped, and % mapped:\t$sampleID\t$n_mapped\t$n_unmapped\t$pct"

fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing output files:"
    [[ "$deinterleave" = true ]] && ls -lh "$R1_mapped" "$R2_mapped" "$R1_unmapped" "$R2_unmapped"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
