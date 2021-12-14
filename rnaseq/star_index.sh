#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=STAR_index
#SBATCH --output=slurm-STAR-index-%j.out


# SETUP ---------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env

## Strict bash settings
set -euo pipefail

## Help function
Help() {
  echo
  echo "## $0: Index a reference genome FASTA file with STAR."
  echo
  echo "## Syntax: $0 -i <input-FASTA> -o <output-dir> [ -a <ref-annotation> ] [ -s <index-size> ] [-sh]"
  echo
  echo "## Options:"
  echo "## -h         Print this help message"
  echo "## -i STR     Input reference FASTA file (REQUIRED)"
  echo "## -o STR     Output dir for index files (REQUIRED)"
  echo "## -a STR     Reference annotation (GFF/GTF) file (default: no GFF/GTF, but this is not recommended)"
  echo "## -s INT     Index size (default: automatically determined from genome size)"
  echo
  echo "## Example: $0 -i refdata/my_genome.fa -o refdata/star_index -a refdata/my_genome.gff"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Hardcoded parameters
READLEN=150

## Option defaults
ref_fa=""
index_dir=""
index_size=""
gff=""

## Parse command-line options
while getopts ':i:a:o:s:h' flag; do
  case "${flag}" in
  i) ref_fa="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  o) index_dir="$OPTARG" ;;
  s) index_size="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Report
echo "## Starting script star_index.sh"
date
echo
echo "## Input FASTA file:             $ref_fa"
[[ "$gff" != "" ]] && echo "## Input GFF file:               $gff"
echo "## Genome index dir (output):    $index_dir"
echo

## Check inputs
[[ ! -f "$ref_fa" ]] && echo "## ERROR: Input FASTA (-i) $ref_fa does not exist" >&2 && exit 1
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && echo "## ERROR: Input file GFF (-a) $gff does not exist" >&2 && exit 1

## Determine index size
if [ "$index_size" = "" ]; then
    echo "## Automatically determining index size..."
    genome_size=$(grep -v "^>" "$ref_fa" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    echo "## Genome size: $genome_size / Index size: $index_size"
else
    echo "## STAR index size:              $index_size"
fi

## If a GFF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then
    gff_arg="--sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent"
else
    gff_arg=""
fi

## Overhang length should be read length minus 1
overhang=$(( READLEN - 1 ))

## Make output dir if needed
mkdir -p "$index_dir"


# RUN STAR ---------------------------------------------------------------------
echo -e "-----------------------\n"

## STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $ref_fa = *gz ]]; then
    echo "## Unzipping gzipped FASTA file..."
    ref_fa_unzip=${ref_fa/.gz/}
    gunzip -c "$ref_fa" > "$ref_fa_unzip"
    ref_fa="$ref_fa_unzip"
fi

echo -e "\n## Indexing genome with STAR...."
STAR --runMode genomeGenerate \
     --genomeDir "$index_dir" \
     --genomeFastaFiles "$ref_fa" \
     --genomeSAindexNbases "$index_size" \
     --sjdbOverhang "$overhang" \
     --runThreadN "$SLURM_CPUS_ON_NODE" ${gff_arg}


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$index_dir"

echo -e "\n## Done with script star_index.sh"
date
