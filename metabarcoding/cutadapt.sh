#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=180
#SBATCH --output=slurm-cutadapt-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

## Help function
Help() {
    echo
    echo "## $0: Run Cutadapt for a single pair of FASTQ files."
    echo
    echo "## Syntax: $0 -i <input-R1-FASTQ> -o <output-dir> ..."
    echo
    echo "## Required options:"
    echo "## -i STR     Input R1 FASTQ file (corresponding R2 will be inferred)"
    echo "## -o STR     Output dir"
    echo
    echo "## Other options:"
    echo "## -f STR     Forward primer sequence"
    echo "## -r STR     Reverse primer sequence"
    echo "## -p STR     File with primer sequences, one pair per line"
    echo "## -d         Don't discard untrimmed sequences (default: discard)"
    echo "## -h         Print this help message"
    echo
    echo "## Example: $0 -i data/sample1_R1.fastq.gz -o results/cutadapt -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "## When you have multiple primer pairs, specify a primer file with -p."
    echo "## The script will compute and use the reverse complements of both primers."
    echo
}

# SETUP ------------------------------------------------------------------------
## Report
echo -e "\n## Starting script cutadapt.sh"
date
echo

## Load software
module load python/3.6-conda5.2  # Load OSC conda module
source activate /users/PAS0471/osu5685/.conda/envs/cutadaptenv # Activate cutadapt environment

## Bash strict settings
set -euo pipefail

## Option defaults
R1_in=""
outdir=""
primer_f=""
primer_r=""
discard_untrimmed=true
primer_file=false

## Parse command-line options
while getopts ':i:o:f:r:p:dh' flag; do
    case "${flag}" in
    i) R1_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    f) primer_f="$OPTARG" ;;
    r) primer_r="$OPTARG" ;;
    d) discard_untrimmed=false ;;
    p) primer_file="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Other parameters
n_cores=$SLURM_CPUS_PER_TASK

## Test input/options
[[ $R1_in = "" ]] && echo "## $0: ERROR: No input FASTQ file (-i) provided" >&2 && exit 1
[[ $outdir = "" ]] && echo "## $0: ERROR: No output dir (-o) provided" >&2 && exit 1

## Report
echo "## Input FASTQ file R1 (-i):   $R1_in"
echo "## Output dir (-o):            $outdir"
echo

## Get primers
primer_arg=""

if [ "$primer_file" = "false" ]; then

    echo "## Using forward and reverse primers provided as arguments..."

    [[ $primer_f = "" ]] && echo "## $0: ERROR: No forward primer (-f) provided" >&2 && exit 1
    [[ $primer_r = "" ]] && echo "## $0: ERROR: No reverse primer (-r) provided" >&2 && exit 1

    primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
    primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

    primer_arg="-a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"

    echo "## Forward primer (-f): $primer_f"
    echo "## Reverse primer (-r): $primer_r"
    echo "## Forward primer - reverse complement: $primer_f_rc"
    echo "## Reverse primer - reverse complement: $primer_r_rc"

else
    echo "## Using primer file $primer_file to read primers..."

    [[ ! -f "$primer_file" ]] && echo -e "\n## $0: ERROR: Primer file $primer_file not found\n" >&2 && exit 1

    while read -r primer_f primer_r; do

        primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
        primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

        echo
        echo "## Forward primer (-f): $primer_f"
        echo "## Reverse primer (-r): $primer_r"
        echo "## Forward primer - reverse complement: $primer_f_rc"
        echo "## Reverse primer - reverse complement: $primer_r_rc"

        primer_arg="$primer_arg -a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"
        primer_arg=$(echo "$primer_arg" | sed -E 's/^ +//') # Remove leading whitespace

    done <"$primer_file"

    echo
fi

## Define cutadapt options
# Including --pair-filter=any is here because cutadapt complains when including empty variable somehow
# "--discard-untrimmed": Remove pairs with no primer found
# "--pair-filter=any": Remove pair if one read is filtered (=Default)
if [ $discard_untrimmed = "true" ]; then
    options="--discard-untrimmed --pair-filter=any"
else
    options="--pair-filter=any"
fi

## Determine input dir and R2 file
indir=$(dirname R1_in)
R2_in=${R1_in/R1/R2}
R1_basename=$(basename "$R1_in")
R2_basename=$(basename "$R2_in")
sample_id=${R1_basename%%_R1}

## Report
echo "## Input FASTQ file R2 (inferred): $R2_in"
echo "## Primer argument: $primer_arg"
echo "## Discard untrimmed (-d): $discard_untrimmed"
echo -e "-----------------------------\n"

## Test input
[[ ! -f "$R1_in" ]] && echo -e "\n## $0: ERROR: Input FASTQ file $R1_in not found\n" >&2 && exit 1
[[ ! -f "$R2_in" ]] && echo -e "\n## $0: ERROR: Input FASTQ file $R2_in not found\n" >&2 && exit 1
[[ "$indir" = "$outdir" ]] && echo "## $0: ERROR: Input dir should not be the same as output dir" >&2 && exit 1

## Create output directory if it doesn't already exist
mkdir -p "$outdir"


# RUN CUTADAPT --------------------------------------------------------
echo "## Running cutadapt..."

cutadapt $primer_arg $options \
    --cores "$n_cores" \
    --output "$outdir"/"$R1_basename" \
    --paired-output "$outdir"/"$R2_basename" \
    "$R1_in" "$R2_in"


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir"/"$sample_id"*

echo -e "\n## Done with script cutadapt.sh."
date
echo
