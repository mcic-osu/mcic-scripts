#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=star_index
#SBATCH --output=slurm-star_index-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "               INDEX A GENOME OR TRANSCRIPTOME WITH STAR"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-FASTA> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--fasta      <file>  Input FASTA file"
    echo "  -o/--outdir     <dir>   Output dir for index files (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --annot         <file>  Reference annotation (GFF/GFF3/GTF) file (GTF preferred)  [default: no annotation file, but this is not recommended]"
    echo "  --index_size    <int>   Index size                                  [default: 'auto' => automatically determined from genome size]"
    echo "  --read_len      <int>   Read length (only applies with --annot)     [default: '150' (bp)]"
    echo "  --overhang      <int>   Overhang (only applies with --annot)        [default: 'auto' => read length minus 1]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to STAR"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for STAR and exit"
    echo "  -v/--version            Print the version of STAR and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i refdata/my_genome.fa -o refdata/star_index -a refdata/my_genome.gff"
    echo
    echo "NOTES:"
    echo "  The script will check how much memory has been allocated to the SLURM job (default: 64GB),"
    echo "  and pass that to STAR via the 'limitGenomeGenerateRAM argument'."
    echo "  When allocating more memory to the SLURM job,"
    echo "  wich can be necessary for large genomes, this will therefore be passed to STAR as well."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - https://github.com/alexdobin/STAR"
    echo "  - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    module load miniconda3/4.12.0-py39
    source activate /fs/project/PAS0471/jelmer/conda/star-2.7.10a  # NOTE: This env includes samtools
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    STAR --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    STAR --help
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
# Option defaults
index_size="auto"
overhang="auto"
read_len=150
mem_bytes=4000000000

debug=false
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
fasta=""
outdir=""
annot=""
more_args=""

# Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --fasta )      shift && fasta=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --annot )           shift && annot=$1 ;;
        --index_size )      shift && index_size=$1 ;;
        --read_len )        shift && read_len=$1 ;;
        --overhang )        shift && overhang=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
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
Load_software
Set_threads

# Determine memory
[[ "$slurm" = true ]] && mem_bytes=$((SLURM_MEM_PER_NODE * 1000000))

# Check input
[[ "$fasta" = "" ]] && Die "Please specify an input file with -i/--fasta" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$fasta" ]] && Die "Input FASTA file $fasta does not exist"
[[ "$annot" != "" ]] && [[ ! -f "$annot" ]] && Die "Input annotation file $annot does not exist"

# Report - part I
echo "=========================================================================="
echo "                    STARTING SCRIPT STAR_INDEX.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input FASTA file:                 $fasta"
echo "Output dir:                       $outdir"
echo "Number of threads/cores:          $threads"
[[ "$annot" != "" ]] && echo "Read length:                      $read_len"
[[ "$annot" != "" ]] && echo "Input annotation file:                   $annot"
[[ $more_args != "" ]] && echo "Other arguments for STAR:         $more_args"
echo
echo "Listing the input file(s):"
ls -lh "$fasta"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
mkdir -pv "$outdir"/logs

# STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $fasta = *gz ]]; then
    fasta_unzip=${fasta/.gz/}
    if [[ ! -f $fasta_unzip ]]; then
        echo "# Unzipping gzipped FASTA file..."
        gunzip -c "$fasta" > "$fasta_unzip"
    else
        echo "# Unzipped version of the FASTA file already exists:"
        ls -lh "$fasta_unzip"
    fi
    fasta="$fasta_unzip"
fi

# Determine index size
if [ "$index_size" = "auto" ]; then
    echo -e "\n# Automatically determining the index size..."
    genome_size=$(grep -v "^>" "$fasta" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    echo "Genome size (autom. determined):  $genome_size"
    echo "Index size (autom. determined):   $index_size"
else
    echo "Index size:                       $index_size"
fi

# If a GFF file is provided, build the appropriate argument for STAR
if [ "$annot" != "" ]; then
    # Overhang length should be read length minus 1 - only if annot is included
    [[ $overhang = "auto" ]] && overhang=$(( read_len - 1 ))
    annot_arg="--sjdbGTFfile $annot --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $overhang"
    echo "Overhang:                     $overhang"
else
    annot_arg=""
fi

# Report
echo "=========================================================================="
# Print reserved resources
[[ "$slurm" = true ]] && Print_resources

# Run STAR
echo -e "\n# Now indexing with STAR...."
Time STAR \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM "$mem_bytes" \
    --genomeDir "$outdir" \
    --genomeFastaFiles "$fasta" \
    --genomeSAindexNbases "$index_size" \
    --runThreadN "$threads" \
    ${annot_arg}


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing files in the output dir:"
ls -lhd "$PWD"/"$outdir"/*
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
