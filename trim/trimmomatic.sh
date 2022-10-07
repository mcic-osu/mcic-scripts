#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-trimmomatic-%j.out


# ARGS AND PARAMETERS ----------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Trim sequences in paired-end FASTQ files using Trimmomatic"
    echo
    echo "Syntax: $0 -i <input-R1-FASTQ> -o <output dir> ..."
    echo
    echo "Required options:"
    echo "  -i STRING     Input R1 (forward reads) sequence file (name of R2 will be inferred)"
    echo "  -o STRING     Output directory (will be created if needed)"
    echo
    echo "Other options:"
    echo "  -a STRING     Adapter file                       [default: 'none']"
    echo "                Possible values: 'NexteraPE-PE.fa', 'TruSeq2-PE.fa', 'TruSeq3-PE.fa'"
    echo "                Or provide the path to your own FASTA file with adapters, e.g. 'adapters.fa' from BBduk"
    echo "                With the default, no adapter trimming is done."
    echo "  -A STRING     Adapter removal parameters         [default: '2:30:10:2:True']"
    echo "  -h            Print this help message and exit"
    echo "  -p STRING     Trimming parameters for Trimmomatic"
    echo "                                                   [default: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36']"
    echo
    echo "Example command:"
    echo "$0 -i data/fastq/A1_R1.fastq.gz -o results/trimmomatic -a metadata/adapters.fa"
    echo
    echo "To submit to the OSC queue, preface with sbatch:"
    echo "sbatch $0 ..."
    echo
    echo "SLURM parameters in script: '--account=PAS0471 --time=3:00:00 --cpus-per-task=12 --output=slurm-trimmomatic-%j.out"
    echo
    echo "Default SLURM parameters can be overridden when submitting the script, e.g.:"
    echo "sbatch -t 15 $0 ...      (override default time reservation of 3 hours, use 15 minutes instead)"
    echo
    echo "Trimmomatic documentation: http://www.usadellab.org/cms/?page=trimmomatic"
    echo
}

## Option defaults
R1_in=""
outdir=""
adapter_file="NA"
trim_param="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
# Same as example on https://github.com/usadellab/Trimmomatic and https://rpubs.com/ednachiang/MetaG_Pipeline
# Alternatively, example of a much stricter mode: "AVGQUAL:28 LEADING:20 TRAILING:20 MINLEN:36"
adapter_param="2:30:10:2:True"
## Same as example on https://github.com/usadellab/Trimmomatic

## Parse options
while getopts ':i:o:a:A:p:h' flag; do
    case "${flag}" in
        i) R1_in="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        a) adapter_file="$OPTARG" ;;
        A) adapter_param="$OPTARG" ;;
        p) trim_param="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo "## ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
        :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/trimmomatic-env

## Bash strict settings
set -euo pipefail

## Other parameters
n_cores="$SLURM_CPUS_ON_NODE"

## Check input
[[ "$R1_in"  = "" ]] && echo "## ERROR: Please provide R1 input FASTQ file with -i flag" && exit 1
[[ "$outdir"  = "" ]] && echo "## ERROR: Please provide output dir with -o flag" && exit 1
[[ ! -f $R1_in ]] && echo "## ERROR: Input file R1_in ($R1_in) does not exist" && exit 1

## Process parameters
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
R1_basename=$(basename "$R1_in" .fastq.gz)
R2_basename=$(basename "$R2_in" .fastq.gz)
sample_ID=${R1_basename/"$R1_suffix"/}

## Check R2
[[ ! -f $R2_in ]] && echo "## ERROR: Input file R2_in ($R2_in) does not exist" && exit 1
[[ "$R1_in" = "$R2_in" ]] && echo "## ERROR: Input R1 and R2 FASTQ files are the same file" >&2 && exit 1

## Define output files
stats_dir="$outdir/log"
discard_dir="$outdir"/discard                          # Dir for discarded sequences
trimstats_file="$stats_dir"/"$sample_ID".trimstats.txt # File with Trimmomatic stdout

R1_out="$outdir"/"$R1_basename".fastq.gz # Output R1 file
R2_out="$outdir"/"$R2_basename".fastq.gz # Output R2 file

R1_discard=$discard_dir/"$R1_basename"_U1.fastq.gz # Output file for discarded R1 reads
R2_discard=$discard_dir/"$R2_basename"_U2.fastq.gz # Output file for discarded R2 reads

## Adapter parameters argument
if [[ $adapter_file = "NA" ]]; then
    adapter_arg=""
else
    adapter_arg=" ILLUMINACLIP:$adapter_file:$adapter_param"
    # As in the example here https://github.com/usadellab/Trimmomatic
fi

## Trimming parameters argument
trim_arg=" $trim_param"

## Report
echo "## Starting script trimmomatic.sh"
date
echo
echo "## Input R1:                     $R1_in"
echo "## Output dir:                   $outdir"
echo "## Trimming argument:            $trim_arg"
echo "## Adapter argument:             $adapter_arg"
echo
echo "## Input R2:                     $R2_in"
echo "## Output R1:                    $R1_out"
echo "## Output R2:                    $R2_out"
echo "## File for discarded R1 seqs:   $R1_discard"
echo "## File for discarded R2 seqs:   $R2_discard"
echo
echo "## Listing input files:"
ls -lh "$R1_in"
ls -lh "$R2_in"
echo -e "-----------------------------\n"


# MAIN -------------------------------------------------------------------------
## Create output dirs
mkdir -p "$discard_dir"   # Create dir for discarded sequences if it doesn't exist
mkdir -p "$stats_dir"     # Create dir for stdout file if it doesn't exist

## Run Trimmomatic
echo -e "## Starting Trimmomatic run..."
trimmomatic PE \
    -threads "$n_cores" \
    "$R1_in" "$R2_in" \
    "$R1_out" "$R1_discard" \
    "$R2_out" "$R2_discard"${adapter_arg}${trim_arg} \
    2>&1 | tee "$trimstats_file"


# WRAP UP ----------------------------------------------------------------------
nreads_raw=$(zcat "$R1_in" | awk '{ s++ } END{ print s/4 }')
nreads_trim=$(zcat "$R1_out" | awk '{ s++ } END{ print s/4 }')
echo -e "\n## Number of raw / trimmed read-pairs: $nreads_raw / $nreads_trim"

echo -e "\n## Listing output files:"
ls -lh "$R1_out"
ls -lh "$R2_out"

echo -e "\n## Done with script trimmomatic.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
