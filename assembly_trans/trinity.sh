#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=170G
#SBATCH --cpus-per-task=42
#SBATCH --job-name=trinity
#SBATCH --output=slurm-trinity-%j.out

#TODO Use trap to copy files from TMPDIR in case of failure

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "============================================================================"
    echo "                            $0"
    echo "                     Run Trinity to assemble a transcriptome"
    echo "============================================================================"
    echo
    echo "USAGE:"
    echo "For de novo assembly, provide FASTQ files witth '--indir' or '--fofn'"
    echo "  sbatch $0 --indir <dir with FASTQ files> -o <output dir> [...]"
    echo "For de novo assembly, provide FASTQ files witth '--indir' or '--fofn'"
    echo "  sbatch $0 --bam <BAM file> -o <output dir> [...]"
    echo "To get help:"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile                 <dir>  Output assembly FASTA file"
    echo "Specify input with one of these options:"
    echo "  A) --indir                  <dir>   Input dir with FASTQ files (for de novo assembly)"
    echo "  B) --fofn                   <file>  File Of File Names (FOFN) with FASTQ files (for de novo assembly)"
    echo "  C) --bam                    <file>  A BAM file (for genome-guided assembly)"              
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --strandedness              <str>   RNAseq library type (Trinity's 'SS_lib_type' argument): 'RF'/'reverse', 'FR'/'forward', or 'unstranded'   [default: 'RF']"
    echo "  --min_contig_length         <int>   Minimum contig length           [default: 200 (= Trinity default)]"
    echo "  --normalize                 <bool>  Whether to normalize reads      [default: 'true']"
    echo "  --normalize_max_read_cov    <int>   Normalize to this coverage      [default: 200 (= Trinity default)]"
    echo "  --genome_guided_max_intron  <int>   Max intron size for genome-guided assembly  [default: 10000]"
    echo "  --workdir                   <dir>   Manually specify a work (initial output) dir [default: auto]"
    echo "  --use_node_tmpdir                   Instead of the OSC scratch temp dir '$PFSDIR', use the node temp dir '$TMPDIR' as the initial output dir"
    echo "  --more_args                 <str>   Quoted string with additional argument(s) to pass to Trinity"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                            Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                             Run the script in debug mode (print all code)"
    echo "  -h                                  Print this help message and exit"
    echo "  --help                              Print the help for Trinity and exit"
    echo "  -v/--version                        Print the version of Trinity and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --indir data/fastq/ -o results/trinity_denovo"
    echo "  sbatch $0 --fofn fofn.txt -o results/trinity_denovo"
    echo "  sbatch $0 --bam results/all_samples.bam -o results/trinity_genomeguided"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/trinityrnaseq/trinityrnaseq/wiki"
    echo
}

# Load the software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/trinity-2.13.2
    set -u
}

# Print version
Print_version() {
    Load_software
    set +e
    Trinity --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    Trinity --show_full_usage_info
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
#                     CONSTANTS AND DEFAULTS
# ==============================================================================
# Defaults
strandedness="RF"
genome_guided=false
genome_guided_max_intron=10000          # Only for genome-guided assembly
normalize=true && normalize_arg=""
min_contig_length=200
normalize_max_read_cov=200
use_node_tmpdir=false                   # Use the scratch dir $PFSDIR instead

debug=false
dryrun=false && e=""
slurm=true

# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
# Placeholder defaults
indir=""
bam=""
fofn=""
outfile=""
outdir_init=""
more_args=""
mem_gb=4        # Will be changed if this is a SLURM job

# Parse command-line options
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outfile )                shift && outfile=$1 ;;
        --workdir )                     shift && outdir_init=$1 ;;
        --indir )                       shift && indir=$1 ;;
        --bam )                         shift && bam=$1 ;;
        --fofn )                        shift && fofn=$1 ;;
        --strandedness)                 shift && strandedness=$1 ;;
        --min_contig_length )           shift && min_contig_length=$1 ;;
        --normalize )                   shift && normalize=$1 ;;
        --normalize_max_read_cov )      shift && normalize_max_read_cov=$1 ;;
        --genome_guided_max_intron )    shift && genome_guided_max_intron=$1 ;;
        --use_node_tmpdir )             use_node_tmpdir=false ;;
        --more_args )                   shift && more_args=$1 ;;
        -v | --version )                Print_version; exit 0 ;;
        -h )                            Print_help; exit 0 ;;
        --help )                        Print_help_program; exit 0;;
        --dryrun )                      dryrun=true && e="echo ";;
        --debug )                       debug=true ;;
        * )                             Die "Invalid option $1" "$all_args" ;;
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

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ "$outfile" = "" ]] && Die "Please specify an output assembly file with -o/--outfile" "$all_args"
[[ "$normalize" != "true" && "$normalize" != "false" ]] && Die "--normalize should be 'true' or 'false', instead it is $normalize"

# Infer output dir and name of gene2transcript map
outdir=$(dirname "$outfile")
outfile_name=(basename "$outfile")
gene2trans="$outdir"/${outfile_name%.*}.gene2trans

# Create comma-delimited list of FASTQ files:
if [[ "$indir" != "" ]]; then
    [[ ! -d "$indir" ]] && Die "Input dir $indir does not exist"
    R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
    R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')
    nfiles=$(find "$indir" -name "*fastq.gz" | wc -l)
    [[ "$nfiles" -eq 0 ]] && Die "No input files found"
elif [[ "$fofn" != "" ]]; then
    [[ ! -f "$fofn" ]] && Die "Input FOFN $fofn does not exist"
    R1_list=$(grep "_R1"  "$fofn" | tr "\n" "," | sed 's/,$/\n/')
    R2_list=$(grep "_R2"  "$fofn" | tr "\n" "," | sed 's/,$/\n/')
    nfiles=$(grep -c "." "$fofn")
    [[ "$nfiles" -eq 0 ]] && Die "No input files found"
elif [[ "$bam" != "" ]]; then
    [[ ! -f "$bam" ]] && Die "Input BAM file $bam does not exist"
    genome_guided=true
else
    Die "Please specify input files with --indir, --fofn, or --bam"
fi

# Library type argument
if [[ "$strandedness" = "unstranded" ]]; then
    strand_arg=""
elif [[ "$strandedness" = "reverse" ]]; then
    strand_arg="--SS_lib_type RF"
elif [[ "$strandedness" = "forward" ]]; then
    strand_arg="--SS_lib_type FR"
else
    strand_arg="--SS_lib_type $strandedness"
fi

# Initial output dir
if [[ "$slurm" = true ]]; then
    if [[ "$outdir_init" = "" ]]; then
        if [[ "$use_node_tmpdir" = true ]]; then
            outdir_init="$TMPDIR"/trinity_out
            mkdir -p "$TMPDIR"/trinity_out
        else
            outdir_init=$PFSDIR/trinity_out
            mkdir -p "$PFSDIR"/trinity_out
        fi
    fi
else
    [[ "$outdir_init" = "" ]] && outdir_init="$outdir"
fi

# Other args
[[ "$normalize" = "false" ]] && normalize_arg="--no_normalize_reads"
[[ "$slurm" = true ]] && mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TRINITY.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
[[ "$indir" != "" ]] && echo "Input dir:                        $indir"
[[ "$fofn" != "" ]] && echo "FOFN:                             $fofn"
[[ "$bam" != "" ]] && echo "Input BAM file:                   $bam"
echo "Output assembly file:             $outfile"
echo "Output gene-to-transcript map:    $gene2trans"
echo "Initial output (work) dir:        $outdir_init"
echo "Strandedness / strand argument:   $strandedness / $strand_arg"
echo "Minimum contig length:            $min_contig_length"
echo "Normalize reads:                  $normalize"
echo "Max cov. to normalize reads to:   $normalize_max_read_cov"
echo "Number of threads/cores:          $threads"
echo "Memory in GB:                     $mem_gb"
[[ $more_args != "" ]] && echo "Other arguments for Trinity:      $more_args"

if [[ "$genome_guided" = false ]]; then
    echo "Run type:                         de novo"
    echo "Number of input files:            $nfiles"
    echo "# Listing the input file(s):"
    [[ $indir != "" ]] && ls -lh "$indir"/*fastq.gz
    [[ "$fofn" != "" ]] && cat "$fofn" | xargs -I{} ls -lh {}
else
    echo "Run type:                         genome-guided"
    echo "Max intron size:                  $genome_guided_max_intron"
    echo "# Listing the input file(s):"
    ls -lh "$bam"
fi

[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
${e}mkdir -pv "$outdir_init" "$outdir"/logs

# Run Trinity
set +e   # Trinity exit codes are weird?
if [[ "$genome_guided" = false ]]; then
    # De novo
    echo -e "\n# Starting de novo Trinity run..."
    ${e}Time Trinity \
        --seqType fq \
        --left "$R1_list" \
        --right "$R2_list" \
        --min_contig_length "$min_contig_length" \
        --normalize_max_read_cov "$normalize_max_read_cov" \
        $normalize_arg \
        $strand_arg \
        --output "$outdir_init" \
        --max_memory "$mem_gb" \
        --CPU "$threads" \
        $more_args
else
    # Genome-guided
    echo -e "\n# Starting genome-guided Trinity run..."
    ${e}Time Trinity \
        --genome_guided_bam "$bam" \
        --genome_guided_max_intron "$genome_guided_max_intron" \
        --min_contig_length "$min_contig_length" \
        --normalize_max_read_cov "$normalize_max_read_cov" \
        $normalize_arg \
        $strand_arg \
        --output "$outdir_init" \
        --max_memory "$mem_gb" \
        --CPU "$threads"
fi
set -e

#? Use '--full_cleanup' option? ":only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta"

# Copy the final assembly files
if [[ "$outdir_init" != "$outdir" ]]; then
    echo -e "\n# Copying assembly files to the final output dir..."
    ${e}find "$outdir_init"/.. \
        -name "*Trinity*fasta*" \
        -not -name "*reads*" \
        -type f \
        -exec cp -vf {} "$outdir" \;
fi

# Move the final assembly files
if [[ "$genome_guided" = false ]]; then
    ${e}mv -v "$outdir"/trinity_out.Trinity.fasta "$outfile"
    ${e}mv -v "$outdir"/trinity_out.Trinity.fasta.gene_trans_map "$gene2trans"
else
    ${e}mv -v "$outdir"/Trinity-GG.fasta "$outfile"
    ${e}mv -v "$outdir"/Trinity-GG.fasta.gene_trans_map "$gene2trans"
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
    ls -lh "$outfile" "$gene2trans"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
