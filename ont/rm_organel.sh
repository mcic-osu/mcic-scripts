#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=rm_organel
#SBATCH --output=slurm-rm_organel-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "       REMOVE MAPPING SEQUENCES FROM LONG-READ FASTQ FILES"
    echo "       (TYPICALLY USED TO REMOVE ORGANELLE-DERIVED SEQUENCES)"
    echo "======================================================================"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--fq_in   <file>     Input FASTQ file"
    echo "  -o/--fq_out  <file>     Output FASTQ file"
    echo "  -r/--ref     <file>     Reference FASTA file"
    echo "  --seqids     <str>      Comma-separated list of sequence IDs (contigs/scaffold)"
    echo "                          from the reference FASTA file that should be removed"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --keep_sam              Keep the intermediate SAM file    [default: remove]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo
}

## Load software
Load_software() {
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do conda deactivate 2>/dev/null; done
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/seqtk
}

## Print args
Print_args() {
    echo -e "\n# Arguments passed to the script:"
    echo "$*"
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo "For help, run this script with the '-h' / '--help' option"
    echo -e "Exiting\n" >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
keep_sam=false

debug=false
dryrun=false

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
fq_in=""
fq_out=""
ref=""
seqids=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --fq_in )          shift && fq_in=$1 ;;
        -o | --fq_out )         shift && fq_out=$1 ;;
        -r | --ref )            shift && ref=$1 ;;
        --seqids )              shift && seqids=$1 ;;
        --keep_sam )            keep_sam=true ;;    
        -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        --dryrun )              dryrun=true ;;
        --debug )               debug=true ;;
        * )                     Print_args "$all_args"; Die "Invalid option $1" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## FASTQ filename parsing
outdir=$(dirname "$fq_out")
file_id=$(basename "$fq_in" .fastq.gz)

## Intermediate files
scaffold_list="$outdir"/intermed/"$file_id"_selected_scaffolds.txt
ref_organel="$outdir"/intermed/"$file_id"_selected_scaffolds.fasta
sam="$outdir"/intermed/$(basename "$ref_organel" .fasta).sam  # Determined by minimap script

## Bash script settings
set -euo pipefail

## Check input
[[ $fq_in = "" ]] && Print_args "$all_args" && Die "Please specify an input FASTQ file with -i/--fq_in"
[[ $fq_out = "" ]] && Print_args "$all_args" && Die "Please specify an output FASTQ file with -o/--fq_out"
[[ $ref = "" ]] && Print_args "$all_args" && Die "Please specify a reference FASTA with -r/--ref"
[[ $seqids = "" ]] && Print_args "$all_args" && Die "Please specify 1 or more sequence IDs with --seqids"
[[ ! -f $fq_in ]] && Die "Input file $fq_in does not exist"
[[ ! -f $ref ]] && Die "Input file $ref does not exist"


## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT RM_ORGANEL.SH"
date
echo "=========================================================================="
echo "Input FASTQ file:             $fq_in"
echo "Output FASTQ file:            $fq_out"
echo "Reference FASTA file:         $ref"
echo "Scaffolds/contigs to target:  $seqids"
echo
echo "Listing the input files:"
ls -lh "$fq_in" "$ref"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN\n"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    
    ## Create the output directory
    mkdir -p "$outdir"/logs "$outdir"/intermed

    ## Run
    echo "$seqids" | tr "," "\n" > "$scaffold_list"

    echo -e "\nNow running seqtk to extract target sequences from reference..."
    seqtk subseq "$ref" "$scaffold_list" > "$ref_organel"

    echo -e "\nNow running the minimap script to map reads to the target sequences..."
    bash mcic-scripts/map/minimap.sh -i "$fq_in" -r "$ref_organel" -o "$outdir"/intermed
    
    echo -e "\nNow running the samtools script to extract unmapped reads..."
    bash mcic-scripts/map/samtools.sh -i "$sam" --command "fastq -n -f 4 -0 $fq_out"

    echo "======================================================================"
    echo -e "\nNow counting sequences in the in- and output files..."
    nseq_in=$(zcat "$fq_in" | awk '{ s++ } END{ print s/4 }')
    nseq_out=$(zcat "$fq_out" | awk '{ s++ } END{ print s/4 }')
    echo "Nr of input / output sequences in FASTQ files: $nseq_in  /  $nseq_out"

    [[ "$keep_sam" = false ]] && rm -v "$sam"

fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
if [[ "$dryrun" = false ]]; then
    echo
    echo "======================================================================"
    echo -e "\nListing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime | grep -Ev "ba|ex"
fi
echo
echo "## Done with script"
date
