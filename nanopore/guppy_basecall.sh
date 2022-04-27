#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm-guppy_basecall-%j.out

# HELP AND COMMAND-LINE OPTIONS ------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Perform FAST5 basecalling with Guppy"
    echo
    echo "Syntax:      $0 [ -i <input-R1-file> / -I <input-R1-dir> ] -o <output-dir> ..."
    echo 
    echo "Required options:"
    echo "  -i STRING          Input FAST5 file (run Guppy on one file)"
    echo "  -I STRING          Input dir with FAST5 files (run Guppy on all files)"
    echo "        NOTE: _Either_ -i or -I has to be be specified"
    echo "  -o STRING          Output directory"
    echo "  -c STRING          Config file                                 [default: ...]"
    echo
    echo "Other options:"
    echo "  -a STRING          Other arguments to pass to Guppy"
    echo "  -b STRING          Barcode set                                 [default: none - no demultiplexing]"
    echo "  -q INTEGER         Minimum Q-score                             [default: 9]"
    echo "  -h                 Print this help message and exit"
    echo
    echo "- Example:"
    echo "      $0 -i data/FAS09095.fast5 -o results/guppy"
    echo 
    echo "- To submit the OSC queue, preface with 'sbatch':"
    echo "      sbatch $0 ..."
    echo 
    echo "- To override one of the SLURM/SBATCH parameters in the script, add these after 'sbatch' - for example:"
    echo "      sbatch --time=60 $0 ..."
    echo
}

## Option defaults
infile=""
indir=""
outdir=""
config=""
barcode_kit=""
more_args=""
min_qscore=9

## Get command-line options
while getopts 'i:I:o:a:b:c:q:h' flag; do
    case "${flag}" in
        i) infile="$OPTARG" ;;
        I) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        c) config="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        b) barcode_kit="$OPTARG" ;;
        q) min_qscore="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
        :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

[[ "$infile" = "" ]] && [[ "$indir" = "" ]] && echo "## ERROR: Please specify either an input file with -i, or an input dir file with -I" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ "$infile" != "" ]] && [[ ! -f "$infile" ]] && echo "## ERROR: input file R1 $infile does note exist" >&2 && exit 1
[[ "$indir" != "" ]] && [[ ! -d "$indir" ]] && echo "## ERROR: input YAML file R1 $indir does note exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software 
GUPPY_DIR=/fs/project/PAS0471/jelmer/software/guppy-6.0.1
GUPPY_BIN="$GUPPY_DIR"/bin/guppy_basecaller

## Bash strict settings
set -euo pipefail

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"/tmp

## Other parameters
N_CORES=$SLURM_CPUS_PER_TASK

## Build input arg
if [[ $infile != "" ]]; then
    indir=$(dirname "$infile")
    infile_base=$(basename "$infile")
    filelist_file="$outdir"/tmp/"$infile_base".filelist.txt
    echo "$infile_base" > "$filelist_file"
    infile_arg="--input_file_list $filelist_file"
fi

## Build barcode arg
if [[ $barcode_kit != "" ]]; then
    barcode_arg="--barcode_kits $barcode_kit --trim_barcodes"
else
    barcode_arg=""
fi

## Report
echo "## Starting script guppy_basecall.sh..."
date
echo
echo "## Input dir with FASTQ files:   $indir"
echo "## Output dir:                   $outdir"
echo "## ONT config name:              $config"
[[ $barcode_kit != "" ]] && echo "## Barcode kit name:             $barcode_kit"
echo "## Minimum qual score to PASS:   $min_qscore"
echo "## Number of cores (threads):    $N_CORES"
echo -e "--------------------\n"


# RUN GUPPY -------------------------------------------------------------------
$GUPPY_BIN \
    --input_path "$indir" ${infile_arg} \
    --save_path "$outdir" \
    --config "$config" \
    --min_qscore "$min_qscore" \
    --compress_fastq \
    --records_per_fastq 0 \
    --cpu_threads_per_caller "$N_CORES" \
    --num_barcode_threads "$N_CORES" \
    --num_callers 1 ${barcode_arg} ${more_args}


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing files in output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script guppy_basecall.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo

## To print config names for flowcell + kit combs:
#$ guppy_basecaller --print_workflows
## For kit LSK109 and flowcell FLO-MIN106, this is dna_r9.4.1_450bps_hac