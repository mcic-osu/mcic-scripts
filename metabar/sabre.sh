#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=sabre
#SBATCH --output=slurm-sabre-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "      RUN SABRE TO DEMULTIPLEX FASTQ FILES USING INLINE BARCODES"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --R1 <R1 FASTQ file> --R2 <R2 FASTQ file> --barcode-file <barcode file> -o <output dir>"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -f/--R1             <file>  Input R1 FASTQ file (can be gzipped)"
    echo "  -r/--R2             <file>  Input R2 FASTQ file (can be gzipped)"
    echo "  -b/--barcode-file   <file>  Input TSV file with 1 line per barcode and 3 columns: barcode, forward output FASTQ, reverse output FASTQ"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo
    echo "OUTPUT:"
    echo "  - Demultiplexed, gzipped R1 and R2 FASTQ files for each sample in the main output dir"
    echo "  - R1 and R2 files with unassigned reads in a directory 'unassigned' within the output dir"
    echo 
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --R1 data/A_R1.fastq.gz --R2 data/A_R2.fastq.gz --barcode_file barcodes.tsv -o results/sabre"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Docs: https://github.com/najoshi/sabre"
    echo "  - Tutorial: https://astrobiomike.github.io/amplicon/demultiplexing"
    echo
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
#                               SETUP
# ==============================================================================
# Parse command-line args
R1=""
R2=""
barcode_file=""
outdir=""

all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -f | --R1 )             shift && R1=$1 ;;
        -r | --R2 )             shift && R2=$1 ;;
        -b | --barcode-file )   shift && barcode_file=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -h | --help )           Print_help; exit 0 ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Load software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/sabre-1.0

# Bash script settings
set -euo pipefail

# Make paths absolute, because we will cd into the outdir
[[ ! "$R1" =~ ^/ ]] && R1="$PWD"/"$R1"
[[ ! "$R2" =~ ^/ ]] && R2="$PWD"/"$R2"
[[ ! "$barcode_file" =~ ^/ ]] && barcode_file="$PWD"/"$barcode_file"
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"

# Define the output files
R1_unassigned_out="$outdir"/unassigned/unassigned_R1.fastq
R2_unassigned_out="$outdir"/unassigned/unassigned_R2.fastq

# Test the input
[[ "$R1" = "" ]] && Die "Please specify an R1 input FASTQ file with --R1" "$all_args"
[[ "$R2" = "" ]] && Die "Please specify an R2 input FASTQ file with --R2" "$all_args"
[[ "$barcode_file" = "" ]] && Die "Please specify a barcode file --barcode-file" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir --outdir" "$all_args"

[[ ! -f "$R1" ]] && Die "Input R1 file $R1 does not exist"
[[ ! -f "$R2" ]] && Die "Input R2 file $R2 does not exist"
[[ ! -f "$barcode_file" ]] && Die "Input barcode file $barcode_file does not exist"


# Report
echo "========================================================================="
echo "                    STARTING SCRIPT SABRE.SH"
date
echo "========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input R1 FASTQ:                   $R1"
echo "Input R2 FASTQ:                   $R2"
echo "Input barcode file:               $barcode_file"
echo "Output dir:                       $outdir"
echo
echo "# Listing the input file(s):"
ls -lh "$R1" "$R2" "$barcode_file"
echo
echo "# Showing the contents of the barcode file:"
cat "$barcode_file"
echo "========================================================================"


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
mkdir -p "$outdir"/logs "$outdir"/unassigned

# Move into the output dir
cd "$outdir" || exit

# Run Sabre
echo -e "\n# Now running Sabre..."
sabre pe \
    --pe-file1 "$R1" \
    --pe-file2 "$R2" \
    --barcode-file "$barcode_file" \
    --unknown-output1 "$R1_unassigned_out" \
    --unknown-output2 "$R2_unassigned_out"

#! -c, --both-barcodes, Optional flag that indicates that both fastq files have barcodes.

# Zip the output files
echo -e "\n# Now zipping the output FASTQ files..."
gzip -v "$outdir"/*fastq


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Listing files in the output dir:"
ls -lhd "$outdir"/*
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
echo
echo "# Done with script"
date
