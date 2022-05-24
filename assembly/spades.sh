#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=slurm-spades-%j.out

# HELP AND COMMAND-LINE OPTIONS ------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Assemble a genome with SPAdes"
    echo
    echo "Syntax:      $0 -i <input-R1-file> -o <output-dir> ..."
    echo 
    echo "Required options:"
    echo "  -i STRING          Input R1 (forward) FASTQ file (the name of the R2 file will be inferred)"
    echo "  -y STRING          YAML file with input files"
    echo "        NOTE: _Either_ -i or -y has to be be specified"
    echo "  -o STRING          Output directory"
    echo
    echo "Other options:"
    echo "  -a STRING          Other arguments to pass to SPAdes"
    echo "  -m STRING          SPAdes run mode                             [default: default SPAdes]"
    echo "                     Possible values: 'isolate', 'meta', 'metaplasmid, 'metaviral', 'plasmid', 'rna', 'rnaviral'"
    echo "  -k STRING          Comma-separated list of kmer sizes          [default: SPAdes default]"
    echo "  -c                 Run in 'careful' mode (small genomes only)  [default: don't run in careful mode]"
    echo "  -C                 Continue an interrupted run                 [default: start anew]"
    echo "  -h                 Print this help message and exit"
    echo
    echo "- Example:"
    echo "      $0 -i data/A1_R1_001.fastq.gz -o results/spades -m meta"
    echo 
    echo "- To submit the OSC queue, preface with 'sbatch':"
    echo "      sbatch $0 ..."
    echo 
    echo "- To override one of the SLURM/SBATCH parameters in the script, add these after 'sbatch' - for example:"
    echo "      sbatch --time=60 $0 ..."
    echo
    echo "- SPAdes documentation: https://github.com/ablab/spades"
}

## Option defaults
R1=""
infile_yaml=""
outdir=""
mode=""
mode_arg=""
kmers=""
careful=false
careful_arg=""
continue=false
more_args=""

## Get command-line options
while getopts 'i:y:o:m:k:a:Cch' flag; do
    case "${flag}" in
    i) R1="$OPTARG" ;;
    y) infile_yaml="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    k) kmers="$OPTARG" ;;
    m) mode="$OPTARG" ;;
    c) careful=true ;;
    C) continue=true ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

[[ "$R1" = "" ]] && [[ "$infile_yaml" = "" ]] && echo "## ERROR: Please specify either an R1 input file with -i, or an input YAML file with -y" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ "$R1" != "" ]] && [[ ! -f "$R1" ]] && echo "## ERROR: input file R1 $R1 does note exist" >&2 && exit 1
[[ "$infile_yaml" != "" ]] && [[ ! -f "$infile_yaml" ]] && echo "## ERROR: input YAML file R1 $infile_yaml does note exist" >&2 && exit 1


# OTHER SETUP ------------------------------------------------------------------
## Report
echo "## Starting script spades.sh..."
date
echo

## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/spades-env

## Bash strict mode
set -euo pipefail

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                            # Retrieve number of cores
mem=$(( (($SLURM_MEM_PER_NODE / 1000)) - 1))             # Convert memory in MB to GB (and subtract 1)

## Build some arguments to pass to SPAdes
[[ "$mode" != "" ]] && mode_arg="--$mode"
[[ "$careful" = true ]] && careful_arg="--careful"

## Input file arg
if [[ "$R1" != "" ]]; then

    ## Infer R2 filename
    R1_suffix=$(echo "$R1" | sed -E 's/.*(_R?[1-2])[_\.][0-9]+\.fastq\.gz/\1/')
    R2_suffix=${R1_suffix/1/2}
    R2=${R1/$R1_suffix/$R2_suffix}
    [[ "$R1" = "$R2" ]] && echo "## ERROR: input file R1 and R2 are the same: $R1" >&2 && exit 1
    [[ ! -f "$R2" ]] && echo "## ERROR: input file R1 $R2 does note exist" >&2 && exit 1

    infile_arg="--pe1-1 $R1 --pe1-2 $R2"
else
    infile_arg="--dataset $infile_yaml"
fi

## Kmer arg
if [[ "$kmers" != "" ]]; then
    kmer_arg="-k $kmers"
else
    kmer_arg=""
fi

## Report
echo "## Command-line args:"
[[ "$R1" != "" ]] && echo "## Input FASTQ file - R1:            $R1"
echo "## Output dir:                       $outdir"
[[ "$kmers" != "" ]] && echo "## kmer sizes:                       $kmers"
echo "## Using 'careful' setting:          $careful"
echo "## Continuing a previous run:        $continue"
[[ "$mode" != "" ]] && echo "## Running mode:                     $mode"
[[ "$more_args" != "" ]] && echo "## Other SPAdes arguments:           $more_args"
echo
echo "## Other variables and settings:"
[[ "$R1" != "" ]] && echo "## Input FASTQ file - R2:            $R2"
echo "## Input file argument:              $infile_arg"
echo "## Number of cores:                  $n_cores"
echo "## Memory in GB:                     $mem"
if [[ $infile_yaml != "" ]]; then
    echo -e "\n Printing contents of input YAML file:\n"
    cat "$infile_yaml"
    echo
fi
echo -e "-----------------------\n\n"

## Create output dir
mkdir -p "$outdir"


# RUN SPADES -------------------------------------------------------------------
if [ "$continue" = false ]; then
    echo "## Now running SPAdes..."
    spades.py $infile_arg \
        -o "$outdir" \
        -t "$n_cores" \
        -m "$mem" \
        ${kmer_arg} ${mode_arg} ${careful_arg} ${more_args}
else
    echo "## Now resuming SPAdes run..."
    spades.py -o "$outdir" --continue
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir"
echo -e "\n## Done with script spades.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
