#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=spades
#SBATCH --output=slurm-spades-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "        Assemble a genome or transcriptome with Spades"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 [ -i <input R1> / -I <indir> / -F <fofn> ] -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outfile    <file>  Output assembly FASTA file"
    echo "Either -i/--R1, -I/--indir, or -f/--fofn has to be be used to specify the input:"
    echo "  --R1            <file>  Input R1 (forward) FASTQ file (the name of the R2 file will be inferred)"
    echo "  --indir         <dir>   Dir with gzipped FASTQ files"
    echo "  --fofn          <file>  File of file names (fofn): one line per input R1 FASTQ file"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --mode          <str>   Spades run mode                             [default: default Spades]"
    echo "                          Possible values: 'isolate', 'meta', 'metaplasmid, 'metaviral', 'plasmid', 'rna', 'rnaviral'"
    echo "  --kmer_size    <str>    Comma-separated list of kmer sizes          [default: 'auto' => Spades default of auto-selecting kmer sizes]"
    echo "  --ss / --strandedness <str>   Strandedness for RNAseq libraries     [default: 'rf' (reverse)]"
    echo "                          Options: 'rf'/'reverse', 'fr'/'forward', or 'unstranded'"
    echo "  --careful               Run in 'careful' mode (small genomes only)  [default: don't run in careful mode]"
    echo "  --continue              Resume an interrupted run                   [default: start anew]"
    echo "  --use_node_tmpdir       Instead of the OSC scratch temp dir '$PFSDIR', use the node temp dir '$TMPDIR'"
    echo "                          This is for Spades' '--tmp-dir' argument."
    echo "                          However, note that sometimes, there isn't enough space on the '$TMPDIR'..."
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to Spades"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Spades and exit"
    echo "  -v/--version            Print the version of Spades and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/A_R1.fastq.gz -o results/spades"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/ablab/spades"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/spades-3.15.5
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    spades.py --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    spades.py --help
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
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - only used in SLURM jobs
# In some cases, spades needs more than the alotted 1TB in $TMPDIR
SCRATCH_TMPDIR="$PFSDIR"
NODE_TMPDIR="$TMPDIR"

# Option defaults
use_node_tmpdir=false
kmer_size="auto" && kmer_arg=""
careful=false && careful_arg=""
continue=false
strandedness="rf"       # Only applies for mode 'rna'
mem=4                   # Only applies when not running a SLURM job
debug="false"
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
R1=""
indir=""
fofn=""
outfile=""
mode=""
mode_arg=""
tmpdir_arg=""
strand_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )             shift && R1=$1 ;;
        -I | --indir )          shift && indir=$1 ;;
        -f | --fofn )           shift && fofn=$1 ;;
        -o | --outfile )        shift && outfile=$1 ;;
        --kmer_size )           shift && kmer_size=$1 ;;
        --ss | --strandedness)  shift && strandedness=$1 ;;
        --mode )                shift && mode=$1 ;;
        --careful )             careful=true ;;
        --continue )            continue=true ;;
        --use_node_tmpdir )     use_node_tmpdir=true ;;  
        --more_args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0;;
        --debug )               debug=true ;;
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

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
Load_software
Set_threads

# Infer outdir
outdir=$(dirname "$outfile")

# Additional variables
if [[ "$slurm" = true ]]; then
    # Memory - convert MB memory to GB (and subtract 1)
    #mem=$(( (SLURM_MEM_PER_NODE / 1000) - 1))
    mem=$(( (SLURM_MEM_PER_NODE / 1000) + 1500)) #! Testing this, following https://github.com/ablab/spades/issues/494
    
    # Determine what to use as the temporary directory
    if [[ "$use_node_tmpdir" = true ]]; then
        TMPDIR_FINAL="$NODE_TMPDIR"
    else
        TMPDIR_FINAL="$SCRATCH_TMPDIR"
    fi
    tmpdir_arg="--tmp-dir=$TMPDIR_FINAL"
fi

# Build some arguments to pass to SPAdes
[[ "$mode" != "" ]] && mode_arg="--$mode"
[[ "$careful" = true ]] && careful_arg="--careful"
[[ "$kmer_size" != "auto" ]] && kmer_arg="-k $kmer_size"

if [[ "$mode" = rna ]]; then
    if [[ "$strandedness" = "unstranded" ]]; then
        strand_arg=""
    elif [[ "$strandedness" = "reverse" ]]; then
        strand_arg="--ss rf"
    elif [[ "$strandedness" = "forward" ]]; then
        strand_arg="--ss fr"
    else
        strand_arg="--ss $strandedness"
    fi
fi

# Input file arg
if [[ "$R1" != "" ]]; then
    # Infer R2 filename
    file_ext=$(basename "$R1" | sed -E 's/.*(.fastq.gz|.fq.gz)/\1/')
    R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2=${R1/$R1_suffix/$R2_suffix}
    infile_arg="--pe1-1 $R1 --pe1-2 $R2"

    [[ "$R1" = "$R2" ]] && Die "Input file R1 and R2 are the same: $R1"
    [[ ! -f "$R2" ]] && Die "Input file R2 ($R2) does not exist"
else
    yaml="$outdir"/input_filenames.yml
    infile_arg="--dataset $yaml"
fi

# Check input
[[ "$R1" = "" ]] && [[ "$indir" = "" ]] && [[ "$fofn" = "" ]] && Die "Please specify either an R1 input file with -i, an input dir with -I, or an input FOFN with -f"
[[ "$outfile" = "" ]] && Die "Please provide an output file with -o/--outfile"
[[ "$R1" != "" ]] && [[ ! -f "$R1" ]] && Die "Input file R1 $R1 does note exist"
[[ "$fofn" != "" ]] && [[ ! -f "$fofn" ]] && Die "Input FOFN $fofn does note exist"

# Report
echo
echo "=========================================================================="
echo "                   STARTING SCRIPT SPADES.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Output assembly file:             $outfile"
echo "Kmer size(s):                     $kmer_size"
echo "Using 'careful' setting:          $careful"
echo "Continuing a previous run:        $continue"
echo "Number of threads/cores:          $threads"
echo "Memory in GB:                     $mem"
[[ "$R1" != "" ]] && echo "Input FASTQ file - R1:            $R1"
[[ "$mode" = rna ]] && echo "Strandedness / strand argument:   $strandedness / $strand_arg"
[[ "$mode" != "" ]] && echo "Running mode:                     $mode"
[[ $more_args != "" ]] && echo "Other arguments for Spades:       $more_args"
[[ $tmpdir_arg != "" ]] && echo "Temp dir argument:                $tmpdir_arg"
[[ "$R1" != "" ]] && echo "Input FASTQ file - R2:            $R2"

# Create the output directory
mkdir -p "$outdir"/logs

# If an input dir was provided, create a FOFN
if [[ "$indir" != "" ]]; then
    echo -e "\nPreparing a FOFN from files in the input dir..."
    fofn="$outdir"/input_filenames.txt
    ls "$indir"/*_R1*q.gz > "$fofn"
fi

# Report part 2
echo "Input file argument:              $infile_arg"
echo "Number of input files:            $(grep -c "." "$fofn")"
echo "Listing the input file(s):"
[[ "$R1" != "" ]] && ls -lh "$R1" "$R2"
[[ "$fofn" != "" ]] && cat "$fofn" | xargs -I{} ls -lh {}
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               PREPARE YAML
# ==============================================================================
if [[ "$fofn" != "" ]]; then
    echo -e "\n # Now creating the input YAML file..."

# Create template file
    cat > "$yaml".tmp1 <<'_EOF'
    [
        {
            orientation: "fr",
            type: "paired-end",
            right reads: [
            R2_reads
            ],
            left reads: [
            R1_reads
            ]
        },
    ]
_EOF

    # Make dirs absolute if needed
    sed "/^\//d" "$fofn" | sed -E "s@^@$PWD/@" > "$yaml".modlist
    sed -n "/^\//p" "$fofn" >> "$yaml".modlist   # Paths that were already absolute

    # Replace placeholder strings with filenames
    R1=$(sed -e 's/^/"/' -e 's/$/"/' "$yaml".modlist | sed 's/$/,/')
    R2=$(sed -e 's/^/"/' -e 's/$/"/' "$yaml".modlist | sed 's/$/,/' | sed 's/_R1/_R2/')
    awk -v r="$R1" '{gsub(/R1_reads/,r)}1' "$yaml".tmp1 > "$yaml".tmp2
    awk -v r="$R2" '{gsub(/R2_reads/,r)}1' "$yaml".tmp2 > "$yaml"

    # Remove temporary files
    rm "$yaml".tmp1 "$yaml".tmp2 "$yaml".modlist

    echo -e "\n # Printing the contents of the input YAML file:"
    cat "$yaml"
    echo "=========================================================================="
fi


# ==============================================================================
#                               RUN
# ==============================================================================
# Run Spades
if [ "$continue" = false ]; then
    echo -e "\n# Running SPAdes..."
    Time spades.py \
        $infile_arg \
        -o "$outdir" \
        --threads "$threads" \
        --memory "$mem" \
        ${kmer_arg} \
        ${mode_arg} \
        ${strand_arg} \
        ${careful_arg} \
        ${tmpdir_arg} \
        ${more_args}
else
    echo -e "\n# Resuming a previous SPAdes run..."
    Time spades.py \
        -o "$outdir" \
        --continue
fi

echo -e "\n# Copying the SPAdes output assembly"
cp -v "$outdir"/transcripts.fasta "$outfile"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing the output assembly file:"
ls -lh "$PWD"/"$outfile"
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
