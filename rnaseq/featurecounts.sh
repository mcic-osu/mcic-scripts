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
    echo "$0: Use featurecounts to create a matrix with per-gene read counts, from a directory of BAM files."
    echo
    echo "Syntax: $0 -i <input-dir> -o <output-dir> -a <gff/gtf-file> ..."
    echo
    echo "Required options:"
    echo "    -i DIR        Input directory with BAM files"
    echo "    -a FILE       Input reference annotation (GFF/GTF) file"
    echo "    -o FILE       Output file with count matrix (e.g. 'counts.txt')"
    echo
    echo "Other options:"
    echo "    -t STRING     Feature type to count                        [default: 'exon']"
    echo "                  (This should correspond to a value in the 3rd column in the GFF/GTF file)"
    echo "    -g STRING     Identifier of the feature type               [default: 'Name' for GFF, 'gene_id' for GTF]"
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
annot_file=""
g_opt=""            # featureCounts default is 'gene_id'
t_opt=exon          # Same default as featureCounts itself

## Parse command-line options
while getopts ':i:o:a:t:g:h' flag; do
    case "${flag}" in
        i) indir="$OPTARG" ;;
        o) outfile="$OPTARG" ;;
        a) annot_file="$OPTARG" ;;
        g) g_opt="$OPTARG" ;;
        t) t_opt="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done

# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/subread-env
MULTIQC_ENV=/fs/project/PAS0471/jelmer/conda/multiqc-1.12

## Strict bash settings
set -euo pipefail

## Check inputs
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-d) $indir does not exist" >&2 && exit 1
[[ ! -f "$annot_file" ]] && echo "## ERROR: Input annotation file (GFF/GTF) (-a) $annot_file does not exist" >&2 && exit 1

## Get outdir parameters
outdir=$(dirname "$outfile")

## Report
echo
echo "## Starting script featurecounts.sh"
date
echo
echo "## BAM input dir (-i):              $indir"
echo "## Output file (-o):                $outfile"
echo "## Annotation (GTF/GFF) file (-a):  $annot_file"
echo "## Feature type (-t):               $t_opt"

## Annotation format
if [[ "$g_opt" = "" && "$annot_file" =~ .*\.gff3? ]]; then
    echo "## Annotation format is GTF, setting aggregation ID to 'Name'"
    g_opt="Name"
elif [[ "$g_opt" = "" && "$annot_file" =~ .*\.gtf ]]; then
    echo "## Annotation format is GTF, setting aggregation ID to 'gene_id'"
    g_opt="gene_id"
else
    echo "## ERROR: Unknown annotation file format" >&2 && exit 1
fi

## Report
echo "## Aggregation ID (-g):             $g_opt"
echo
echo "## Number of BAM files:             $(find "$indir"/*bam | wc -l)"
echo -e "-------------------\n"


# MAIN -------------------------------------------------------------------------
## Make output dir if needed
mkdir -p "$outdir"

## Run featurecounts
featureCounts \
    -s 2 \
    -p \
    -B \
    -C \
    -F GTF \
    -t "$t_opt" \
    -g "$g_opt" \
    -a "$annot_file" \
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

## Run MultiQC
conda activate "$MULTIQC_ENV"
multiqc --force --interactive "$outdir" -o "$outdir"


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
