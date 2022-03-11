#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=featurecounts
#SBATCH --out=slurm-featurecounts-%j.out


# SETUP ---------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Create a matrix with per-gene read counts for a directory of BAM files."
  echo
  echo "## Syntax: $0 -i <input-FASTA> -o <output-dir> -a <gff-file> ..."
  echo
  echo "## Required options:"
  echo "## -i STR     Input directory with BAM files"
  echo "## -o STR     Output directory"
  echo "## -a STR     Reference annotation (GFF/GTF) file"
  echo "## Other options:"
  echo "## -t STR     Feature type in GFF file to count (default: 'gene')"
  echo "## -g STR     Name of the feature type in the GFF file (default: 'Name')"
  echo "## -h         Print this help message"
  echo
  echo "## Example: $0 -i results/bam -o results/featurecounts -a refdata/my_genome.gff"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/subread-env

## Strict bash settings
set -euo pipefail

## Option defaults
indir=""
outdir=""
gff=""
t_opt=gene
g_opt=Name

## Parse command-line options
while getopts ':i:o:a:t:g:h' flag; do
  case "${flag}" in
  i) indir="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  t) g_opt="$OPTARG" ;;
  g) t_opt="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Process parameters
outfile=$outdir/counts.txt

## Report
echo
echo "## Starting script featurecounts.sh"
date
echo

## Check inputs
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-d) $indir does not exist" >&2 && exit 1
[[ ! -f "$gff" ]] && echo "## ERROR: Input file GFF (-a) $gff does not exist" >&2 && exit 1

## Report
echo "## BAM input dir (-i):              $indir"
echo "## Output dir (-o):                 $outdir"
echo "## Annotation (GTF/GFF) file (-a):  $gff"
echo "## Feature type (-t):               $t_opt"
echo "## Aggregation ID (-g):             $g_opt"
echo
echo "## Output file:                     $outfile"
echo "## Number of BAM files:             $(find "$indir"/*bam | wc -l)"
echo -e "-------------------\n"

## Make output dir if needed
mkdir -p "$outdir"


# RUN FEATURECOUNTS ------------------------------------------------------------
featureCounts \
    -s 2 \
    -p \
    -B \
    -C \
    -t "$t_opt" \
    -g "$g_opt" \
    -a "$gff" \
    -o "$outfile" \
    -T "$SLURM_CPUS_ON_NODE" \
    "$indir"/*bam

## Options used:
#? -s 2  => Reverse-stranded library like TruSeq
#? -p    => Count fragments, not reads (paired-end)
#? -B    => Require both members of a read pair to be aligned
#? -C    => Don't count pairs with discordant mates

## Other possible options:
#? -O    => Assign reads that overlap multiple features
#? -M    => Include multi-mapping reads
#? --minOverlap => Min nr of overlapping bases required for read assignnment (default: 1)


# WRAP UP ----------------------------------------------------------------------
echo -e "\n------------------------\n## Listing output file:"
ls -lh "$outfile"

echo -e "\n## Done with script featurecounts.sh"
date
