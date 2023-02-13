#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=quast
#SBATCH --output=slurm-quast-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "   Run QUAST to check the quality of one or more genome assemblies"
    echo "======================================================================"
    echo "USAGE:"
    echo "  sbatch $0 -o <output-dir> [...] [ --assembly <assembly> | --assembly_dir <input-dir> | assembly1 assembly2 ... ]"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outdir         <dir>   Output dir"
    echo "  To specify the input assembly/assemblies, use one of the following options:"
    echo "    A) -i/--assembly  <file>  Input assembly FASTA file"
    echo "    B) --assembly_dir <dir>   Input dir with assembly FASTA files (extension '.fasta')"
    echo "    C) Pass assembly FASTA files(s) as positional arguments at the end of the command."
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --ref_fa            <file>  Reference genome nucleotide FASTA file      [default: no reference genome]"
    echo "  --ref_annot         <file>  Reference genome annotation (GFF/GTF) file  [default: no reference genome]"
    echo "  --R1                <file>  FASTQ file with forward (R1) Illumina reads [default: no reads]"
    echo "                              (The R2 filename will be inferred)"
    echo "  --reads             <file>  FASTQ file with single-end Illumina reads"
    echo "  --fragmented                QUAST's '--fragmented' option, use for fragmented assemblies"
    echo "  --large                     QUAST's '--large' option, recommended for genomes >100 Mbp"
    echo "  --kmer_stats                QUAST's '--k-mer-stats' option, recommended for genomes >100 Mbp"
    echo "  --run_busco                 Use this flag to run BUSCO within QUAST"
    echo "                              BUSCO is not run by default, because at least for eukaryotes it fails to download the database"
    echo "  --more_args         <str>   Quoted string with additional argument(s) to pass to QUAST"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for Quast and exit"
    echo "  -v/--version                Print the version of Quast and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --assembly results/assembly/my.fasta -o results/quast"
    echo "  sbatch $0 --assembly_dir results/assemblies -o results/quast"
    echo "  sbatch $0 -o results/quast results/racon/assembly.fasta results/smartdenovo/assembly.fasta"
    echo "  sbatch $0 --assembly_dir results/assemblies -o results/quast --more_args '--eukaryote'"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - GitHub repo: https://github.com/ablab/quast"
    echo "  - Manual: https://quast.sourceforge.net/docs/manual.html"
    echo "  - Paper: https://academic.oup.com/bioinformatics/article/29/8/1072/228832"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/quast-5.0.2
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    quast.py --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    quast.py --help
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
kmer_stats=false && kmer_stats_arg=""
is_large=false && is_large_arg=""
is_fragmented=false && is_fragmented_arg=""
run_busco=false        # Not run by default -- at least for eukaryotes, it fails to download the database file

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
assembly=""
assembly_dir=""
declare -a assemblies
outdir=""
ref_fa="" && ref_fa_arg=""
ref_annot="" && ref_annot_arg=""
R1=""
R2=""
reads="" && illumina_reads_arg=""
nanopore="" && nanopore_arg=""

busco_arg=""
more_args=""

# Parse command-line args
all_args="$*"
count=0
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )     shift && outdir=$1 ;;
        -i | --assembly )   shift && assembly=$1 ;;
        --assembly_dir )    shift && assembly_dir=$1 ;;
        --ref_fa )          shift && ref_fa=$1 ;;
        --ref_annot )       shift && ref_annot=$1 ;;
        --R1 )              shift && R1=$1 ;;
        --reads )           shift && reads=$1 ;;
        --nanopore )        shift && nanopore=$1 ;;
        --fragmented )      is_fragmented=true ;;
        --large )           is_large=true ;;
        --kmer_stats )      kmer_stats=true ;;
        --run_busco )       run_busco=true ;;
        --more_args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 assemblies[count]=$1 && count=$(( count + 1 )) ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir"
[[ "$assembly" != "" && ! -f "$assembly" ]] && Die "Input file (--assembly) $assembly does not exist or is not a regular file"
[[ "$assembly_dir" != "" && ! -d "$assembly_dir" ]] && Die "Input dir (--assembly_dir) $assembly_dir does not exist or is not a directory"

# Build argument for assembly input
if [[ $assembly_dir != "" ]]; then
    mapfile -t assemblies < <(find "$assembly_dir" -name "*.fasta")
elif [[ $assembly != "" ]]; then
    assemblies=("$assembly")
elif [[ ${#assemblies[@]} -eq 0 ]]; then
    Die "Please specify input with --assembly_dir, --assembly or positional arguments"
fi

# If the input is a single file, make a separate output dir
if [[ $assembly != "" ]]; then
    sampleID=$(basename "$assembly" | sed -E 's/.fn?as?t?a?//')
    outdir=$outdir/"$sampleID"
fi

# Build input reads arg(s)
if [[ $R1 != "" ]]; then
    file_ext=$(basename "$R1" | sed -E 's/.*(.fastq.gz|.fq.gz)$/\1/')
    R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2=${R1/$R1_suffix/$R2_suffix}
    illumina_reads_arg="-1 $R1 -2 $R2"
    [[ ! -f "$R1" ]] && Die "Input R1 reads file $R1 does not exist"
    [[ ! -f "$R2" ]] && Die "Input R2 reads file $R2 does not exist"
    [[ "$R1" = "$R2" ]] && Die "R1 and R2 FASTQ files are the same file: $R1"
elif [[ $reads != "" ]]; then
    illumina_reads_arg="--single $reads"
    [[ ! -f "$reads" ]] && Die "Input single-end reads file $reads does not exist"
fi

if [[ $nanopore != "" ]]; then
    nanopore_arg="--nanopore $nanopore"
    [[ ! -f "$nanopore" ]] && Die "Input Nanopore reads file $nanopore does not exist"
fi

# Build other arguments
[[ "$ref_fa" != "" ]] && ref_fa_arg="-r $ref_fa"
[[ "$ref_annot" != "" ]] && ref_annot_arg="--features $ref_annot"
[[ "$is_fragmented" = true ]] && is_fragmented_arg="--fragmented"
[[ "$is_large" = true ]] && is_large_arg="--large"
[[ "$kmer_stats" = true ]] && kmer_stats_arg="--k-mer-stats"
[[ "$run_busco" = true ]] && busco_arg="--conserved-genes-finding"

# Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT QUAST.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Output dir:                           $outdir"
echo "Assemblies:                           ${assemblies[*]}"
echo "Number of assemblies:                 ${#assemblies[@]}"
[[ $assembly != "" ]] && echo "Input assembly FASTA:                 $assembly"
[[ $assembly_dir != "" ]] && echo "Input dir with assemblies:            $assembly_dir"
[[ $ref_fa != "" ]] && echo "Reference FASTA file:                 $ref_fa"
[[ $ref_annot != "" ]] && echo "Reference annotation file:            $ref_annot"
[[ $R1 != "" ]] && echo "R1 FASTQ file:                        $R1"
[[ $R2 != "" ]] && echo "R2 FASTQ file:                        $R2"
[[ $reads != "" ]] && echo "Single-end FASTQ file:                $reads"
[[ $nanopore != "" ]] && echo "Nanopore FASTQ file:                  $nanopore"
[[ $more_args != "" ]] && echo "Other arguments to pass to QUAST:     $more_args"
echo "Fragmented assembly:                  $is_fragmented"
echo "Large assembly (>100 Mbp):            $is_large"
echo "Use kmer stats:                       $kmer_stats"
echo -e "\n# Listing the input files:"
[[ "$ref_fa" != "" ]] && ls -lh "$ref_fa"
[[ "$ref_annot" != "" ]] && ls -lh "$ref_annot"
[[ "$R1" != "" ]] && ls -lh "$R1" "$R2"
[[ "$reads" != "" ]] && ls -lh "$reads" 
[[ "$nanopore" != "" ]] && ls -lh "$nanopore"
ls -lh "${assemblies[@]}"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Now creating the output directories..."
${e}mkdir -pv "$outdir"/logs

# Run
echo -e "\n# Now running QUAST..."
${e}Time quast.py \
    --threads "$threads" \
    --output-dir "$outdir" \
    --glimmer \
    --circos \
    --no-snps \
    $busco_arg \
    $ref_fa_arg \
    $ref_annot_arg \
    $illumina_reads_arg \
    $nanopore_arg \
    $is_fragmented_arg \
    $is_large_arg \
    $kmer_stats_arg \
    $more_args \
    "${assemblies[@]}"

#? - glimmer: Acticate gene-finding, use Glimmer instead of GeneMark tools (trouble running those due to licensing issues)
#? - conserved-genes-finding => runs Busco

#? -k  --k-mer-stats   Compute k-mer-based quality metrics (recommended for large genomes)
#?                     This may significantly increase memory and time consumption on large genomes


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
    [[ "$slurm" = true ]] && Resource_usage
fi
echo
echo "# Done with script"
date
