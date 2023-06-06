#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=kraken
#SBATCH --output=slurm-kraken-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "===================================================================================="
    echo "                  $0"
    echo "Run Kraken2 to assign taxonomy to sequences in a FASTA/FASTQ/pair of FASTQ file(s)"
    echo "===================================================================================="
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-file> -o <output-dir> -d <kraken-db-dir> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile             <file>  Input sequence file (FASTA, single-end FASTQ, or R1 from paired-end FASTQ)"
    echo "                                    - If an R1 paired-end FASTQ file is provided, the name of the R2 file will be inferred"
    echo "                                    - FASTA files should be unzipped; FASTQ files should be gzipped"
    echo "  -o/--outdir             <dir>   Output directory"
    echo "  -d/--db-dir             <dir>   Directory with an existing Kraken database"
    echo "                                    - A few databases are available at: /fs/ess/PAS0471/jelmer/refdata/kraken"
    echo "                                    - Kraken databases can be downloaded from: https://benlangmead.github.io/aws-indexes/k2"
    echo "                                    - Finally, you can use 'kraken-build-custom-db.sh' or 'kraken-build-std-db.sh' to create a database"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --confidence            <num>   Confidence required for assignment: number between 0 and 1            [default: 0.5]"
    echo "  --minimum-base-quality  <int>   Base quality Phred score required for use of a base in assignment     [default: 0]"
    echo "                                    - NOTE: When not 0, output sequences contain 'x's for masked bases"
    echo "  --memory-mapping                Don't load the full database into RAM memory                          [default: load into memory]"
    echo "                                    - Considerably lower, but useful/needed with very large databases"
    echo "  --use-names                     Add taxonomic names to the Kraken 'main' output file                  [default: don't add]"
    echo "                                    - NOTE: This option is not compatible with Krona plotting"
    echo "  --single-end                    FASTQ files are single-end                                            [default: paired-end]"
    echo "  --classified-out                Write 'classified' sequences to file in '<outdir>/classified' dir     [default: don't write]"
    echo "  --unclassified-out              Write 'unclassified' sequences to file in '<outdir>/unclassified' dir [default: don't write]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Kraken and exit"
    echo "  -v/--version            Print the version of Kraken and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/A1_R1.fastq.gz -o results/kraken -d /fs/project/PAS0471/jelmer/refdata/kraken/std"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/kraken2-2.1.2
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    kraken2 --version
    set -e
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
min_conf=0.5
min_q=0
add_names=false && names_arg=""
use_ram=true && mem_map_arg=""
single_end=false

slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder variables
infile=""
outdir=""
db_dir=""
write_classif="" && class_out_arg=""
write_unclassif="" && unclass_out_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )             shift && infile=$1 ;;
        -o | --outdir )             shift && outdir=$1 ;;
        -d | --db-dir )             shift && db_dir=$1 ;;
        -c | --confidence )         shift && min_conf=$1 ;;
        -q | --minimum-base-quality )   shift && min_q=$1 ;;
        --use-names )               add_names=true ;;
        -s | --single-end )         single_end=true ;;
        -w | --classified-out )     write_classif=true ;;
        -W | --unclassified-out )   write_unclassif=true ;;   
        --memory-mapping )          use_ram=false ;;
        --more-args )               shift && more_args=$1 ;;
        -v | --version )            Print_version; exit 0 ;;
        -h )                        Print_help; exit 0 ;;
        --help )                    Print_help_program; exit 0;;
        * )                         Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
Load_software
Set_threads

# Check input
[[ "$infile" = "" ]] && Die "Must specify input file with -i" "$all_args" 
[[ "$db_dir" = "" ]] && Die "Must specify a Kraken DB dir with -d" "$all_args"
[[ "$outdir" = "" ]] && Die "Must specify an output dir with -o" "$all_args"
[[ "$min_conf" = "" ]] && Die "Min confidence is not set" "$all_args"
[[ "$min_q" = "" ]] && Die "Min quality is not set" "$all_args"
[[ ! -f "$infile" ]] && Die "Input file $infile does note exist"
[[ ! -d "$db_dir" ]] && Die "Kraken DB dir $db_dir does note exist"

# Report - part 1
echo
echo "=========================================================================="
echo "                     STARTING SCRIPT KRAKEN.SH"
date
echo "=========================================================================="
echo "Input file:                     $infile"
echo "Output dir:                     $outdir"
echo "Kraken db dir:                  $db_dir"
echo
echo "Add tax. names:                 $add_names"
echo "Min. base qual:                 $min_q"
echo "Min. confidence:                $min_conf"
echo "Use RAM to load database:       $use_ram"
echo

# Build Kraken args (leave space after!)
# RAM
[[ "$use_ram" = false ]] && mem_map_arg="--memory-mapping "
# Add tax. names or not -- when adding names, can't use the output for Krona plotting  
[[ "$add_names" = true ]] && names_arg="--use-names "

# Make sure input file argument is correct based on file type 
if [[ "$infile" =~ \.fa?s?t?q.gz$ ]]; then

    R1_in="$infile"
    R1_basename=$(basename "$R1_in" | sed -E 's/\.fa?s?t?q\.gz//')
    R1_suffix=$(echo "$R1_basename" | sed -E 's/.*(_R?[12]).*/\1/')
    
    if [[ "$single_end" = false ]]; then
        echo "Input type is:                paired-end FASTQ"
        file_type=pe
        R2_suffix=${R1_suffix/1/2}
        R2_in=${R1_in/$R1_suffix/$R2_suffix}
        sample_ID=${R1_basename/"$R1_suffix"/}
        infile_arg="--gzip-compressed --paired $R1_in $R2_in"

        echo "Input FASTQ file - R1:        $R1_in"
        echo "Input FASTQ file - R2:        $R2_in"

        [[ ! -f "$R2_in" ]] && Die "R2 file $R2_in does not exist"
        [[ "$R1_in" = "$R2_in" ]] && Die "R1 file $R1_in is the same as R2 file $R2_in"

        if [[ "$write_classif" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID#.fastq "
        fi
        if [[ "$write_unclassif" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID#.fastq "
        fi

    else
        echo "Input type is:                single-end FASTQ"
        file_type=se
        sample_ID=$(basename "$R1_in" .fastq.gz)
        infile_arg="--gzip-compressed $R1_in"

        if [[ "$write_classif" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID.fastq "
        fi
        if [[ "$write_unclassif" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fastq "
        fi
    
    fi

elif [[ "$infile" =~ \.fn?a?s?t?a$ ]]; then
    echo -e "Input type is:              FASTA"
    file_type=fasta
    infile_basename=$(basename "$infile")
    sample_ID=${infile_basename%%.*}
    infile_arg="$infile"

    if [[ "$write_classif" = true ]]; then
        class_out_arg="--classified-out $outdir/classified/$sample_ID.fa "
    fi
    if [[ "$write_unclassif" = true ]]; then
        unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fa "
    fi

else
    Die "Unknown input file type"
fi

# Define output text files
outfile_main="$outdir"/"$sample_ID".main.txt
outfile_report="$outdir"/"$sample_ID".report.txt

# Report
echo "Number of threads/cores:        $threads"
echo "Sample ID (inferred):           $sample_ID"
echo "Output file - main:             $outfile_main"
echo "Output file - report:           $outfile_report"
[[ "$write_classif" = true ]] && echo "Writing classified sequences:   $class_out_arg"
[[ "$write_unclassif" = true ]] && echo "Writing unclassified sequences: $unclass_out_arg"
echo
echo "Listing the input file(s):"
ls -lh "$infile"
echo "=========================================================================="
echo

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create output dirs
echo -e "\n# Creating the output directories..."
[[ "$write_classif" = true ]] && mkdir -pv "$outdir"/classified
[[ "$write_unclassif" = true ]] && mkdir -pv "$outdir"/unclassified
mkdir -pv "$outdir"/logs

# Run Kraken
echo -e "\n# Starting Kraken2 run..."
Time kraken2 ${more_args}${names_arg}--threads "$threads" \
    ${mem_map_arg}--minimum-base-quality "$min_q" \
    --confidence "$min_conf" \
    --report-minimizer-data \
    ${unclass_out_arg}--db "$db_dir" \
    ${class_out_arg}--report "$outfile_report" \
    ${infile_arg}> "$outfile_main"

#? report-minimizer-data: see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information

# Rename and zip FASTQ files - only implemented for PE FASTQ
if  [[ "$file_type" = "pe" ]]; then
    if [[ "$write_classif" = true ]]; then
        echo -e "\n# Zipping FASTQ files with classified reads..."
        mv "$outdir"/classified/"$sample_ID"_1.fastq "$outdir"/classified/"$sample_ID"_R1.fastq
        gzip -f "$outdir"/classified/"$sample_ID"_R1.fastq

        mv "$outdir"/classified/"$sample_ID"_2.fastq "$outdir"/classified/"$sample_ID"_R2.fastq
        gzip -f "$outdir"/classified/"$sample_ID"_R2.fastq
    fi

    if [[ "$write_unclassif" = true ]]; then
        echo -e "\n# Zipping FASTQ files with unclassified reads..."
        mv "$outdir"/unclassified/"$sample_ID"_1.fastq "$outdir"/unclassified/"$sample_ID"_R1.fastq
        gzip -f "$outdir"/unclassified/"$sample_ID"_R1.fastq

        mv "$outdir"/unclassified/"$sample_ID"_2.fastq "$outdir"/unclassified/"$sample_ID"_R2.fastq
        gzip -f "$outdir"/unclassified/"$sample_ID"_R2.fastq
    fi
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "\n# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing output files:"
ls -lh "$outfile_main" "$outfile_report"
[[ "$write_classif" = true ]] && ls -lh "$outdir"/classified/"$sample_ID"*
[[ "$write_unclassif" = true ]] && ls -lh "$outdir"/unclassified/"$sample_ID"*
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
