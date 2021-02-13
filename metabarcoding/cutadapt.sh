#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=60
#SBATCH --account=PAS0471

set -e -u -o pipefail # Run bash in "safe mode", basically

# SETUP --------------------------------------------------------

# Software at OSC:
module load python/3.6-conda5.2                                # Load conda module
source activate /users/PAS0471/osu5685/.conda/envs/cutadaptenv # Activate cutadapt environment

# Help:
Help() {
  # Display Help
  echo
  echo "## $0: Run cutadapt on all .fastq.gz files found recursively under a specified input dir."
  echo
  echo "## Syntax: $0 -i <input-dir> -o <output-dir> -f <forward-primer> -r <reverse-primer> [-d] [-h]"
  echo "## Options:"
  echo "## -h     Print help."
  echo "## -i     Input dir (REQUIRED)"
  echo "## -o     Output dir (REQUIRED)"
  echo "## -f     Forward primer (REQUIRED)"
  echo "## -r     Reverse primer (REQUIRED)"
  echo "## -d     Don't discard untrimmed (default: discard)"
  echo "## Example: $0 -i data/fastq/raw -o data/fastq/trimmed -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT [-h]"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo "## The script will compute and use the reverse complements of both primers."
  echo
}

# Option defaults:
indir=""
outdir=""
primer_f=""
primer_r=""
discard_untrimmed=true
primer_file=false

# Get command-line options:
while getopts ':i:o:f:r:p:dh' flag; do
  case "${flag}" in
  i) indir="$OPTARG" ;;
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

[[ $indir = "" ]] && echo "## $0: ERROR: No input dir (-i) provided" >&2 && exit 1
[[ $outdir = "" ]] && echo "## $0: ERROR: No output dir (-o) provided" >&2 && exit 1

# Report:
echo -e "\n## Starting cutadapt script."
date
echo
echo "## Using the following parameters:"
echo "## Input dir (-i): $indir"
echo "## Output dir (-o): $outdir"

# Get primers:
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
  
  done < "$primer_file"

  echo
fi

# Set options:
# Including --pair-filter=any is here because cutadapt complains when including empty variable somehow
# "--discard-untrimmed": Remove pairs with no primer found
# "--pair-filter=any": Remove pair if one read is filtered (=Default)
if [ $discard_untrimmed = "true" ]; then
  options="--discard-untrimmed --pair-filter=any"
else
  options="--pair-filter=any"
fi

# Report:
echo "## Primer argument: $primer_arg"
echo "## Discard untrimmed (-d): $discard_untrimmed"

# Test:
[[ ! -d "$indir" ]] && echo -e "\n## $0: ERROR: Input directory not found\n" >&2 && exit 1
[[ $(find . -name "*fastq.gz" | wc -c) = 0 ]] && echo -e "\n## $0: ERROR: No fastq files found\n" >&2 && exit 1

# Create output directory if it doesn't already exist:
mkdir -p "$outdir"

# RUN CUTADAPT --------------------------------------------------------

echo -e "\n## Looping through input files...\n"

shopt -s globstar nullglob # Turn on recursive globbing

for R1 in "$indir"/**/*_R1*.fastq.gz; do

  R2=${R1/_R1_/_R2_}

  R1_basename=$(basename "$R1")
  R2_basename=$(basename "$R2")

  # Report input files:
  echo -e "\n------------------------------------------------------------\n"
  echo "## R1 input file:"
  ls -lh "$R1"
  echo "## R2 input file:"
  ls -lh "$R2"

  # Trim:
  echo -e "\n\n## Running cutadapt..."

  cutadapt $primer_arg $options \
    -o "$outdir"/"$R1_basename" -p "$outdir"/"$R2_basename" "$R1" "$R2"

  # Options:
  # "-a"/"-A": Primers for R1/R2
done

# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir"

echo -e "\n## Done with cutadapt script."
date
