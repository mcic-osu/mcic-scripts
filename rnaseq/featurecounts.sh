#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=featurecounts
#SBATCH --out=slurm-featurecounts-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Create a matrix with per-gene read counts for a directory of BAM files."
  echo
  echo "Syntax: $0 -i <input-FASTA> -o <output-dir> -a <gff-file> ..."
  echo
  echo "Required options:"
  echo "    -i DIR        Input directory with BAM files"
  echo "    -a FILE       Input reference annotation (GFF/GTF) file"
  echo "    -o FILE       Output file with count matrix (e.g. 'counts.txt')"
  echo
  echo "Other options:"
  echo "    -g STRING     Feature type to count                        [default: 'gene']"
  echo "                  (This should correspond to a value in the 3rd column in the GFF/GTF file)"
  echo "    -t STRING     Identifier of the feature type               [default: 'Name']"
  echo "                  (This should correspond to the key for the desired feature type (e.g. gene) in the last column in the GFF/GTF file)"
  echo "    -h            Print this help message and exit"
  echo
  echo "Example: $        0 -i results/bam -o results/featurecounts -a refdata/my_genome.gff"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
indir=""
outfile=""
gff=""
t_opt=gene
g_opt=Name

## Parse command-line options
while getopts ':i:o:a:t:g:h' flag; do
  case "${flag}" in
  i) indir="$OPTARG" ;;
  o) outfile="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  g) g_opt="$OPTARG" ;;
  t) t_opt="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
  :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

# SETUP ---------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/subread-env

## Strict bash settings
set -euo pipefail

## Process parameters
outdir=$(dirname "$outfile")

## Check inputs
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-d) $indir does not exist" >&2 && exit 1
[[ ! -f "$gff" ]] && echo "## ERROR: Input file GFF (-a) $gff does not exist" >&2 && exit 1

## Report
echo
echo "## Starting script featurecounts.sh"
date
echo
echo "## BAM input dir (-i):              $indir"
echo "## Output file (-o):                $outfile"
echo "## Annotation (GTF/GFF) file (-a):  $gff"
echo "## Feature type (-t):               $t_opt"
echo "## Aggregation ID (-g):             $g_opt"
echo
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
    -F GTF \
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
echo -e "\n------------------------"
echo "## Listing output file:"
ls -lh "$outfile"
echo
echo "## Done with script featurecounts.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
