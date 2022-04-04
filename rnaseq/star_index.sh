#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=STAR_index
#SBATCH --output=slurm-STAR-index-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Index a reference genome FASTA file with STAR."
  echo
  echo "## Syntax: $0 -i <input-FASTA> -o <output-dir> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING     Input reference FASTA file"
  echo "## -o STRING     Output dir for index files"
  echo
  echo "## Other options:"
  echo "## -a STRING      Reference annotation (GFF/GTF) file [default: no GFF/GTF, but this is not recommended]"
  echo "## -s INTEGER     Index size                          [default: automatically determined from genome size]"
  echo "## -r INTEGER     Read length                         [default: '150' (bp)]"
  echo "## -v INTEGER     Overhang                            [default: read length minus 1]"
  echo "## -h             Print this help message and exit"
  echo
  echo "## Example: $0 -i refdata/my_genome.fa -o refdata/star_index -a refdata/my_genome.gff"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
ref_fa=""
index_dir=""
index_size=""
gff=""
overhang=""
read_len=150

## Parse command-line options
while getopts ':i:a:o:s:r:v:h' flag; do
  case "${flag}" in
  i) ref_fa="$OPTARG" ;;
  a) gff="$OPTARG" ;;
  o) index_dir="$OPTARG" ;;
  s) index_size="$OPTARG" ;;
  r) read_len="$OPTARG" ;;
  v) overhang="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env

## Strict bash settings
set -euo pipefail

## Report
echo -e "\n## Starting script star_index.sh"
date
echo

## Check inputs
[[ ! -f "$ref_fa" ]] && echo "## ERROR: Input FASTA (-i) $ref_fa does not exist" >&2 && exit 1
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && echo "## ERROR: Input file GFF (-a) $gff does not exist" >&2 && exit 1

## Determine index size
if [ "$index_size" = "" ]; then
    echo "## Automatically determining index size..."
    genome_size=$(grep -v "^>" "$ref_fa" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    echo "## Genome size:                  $genome_size"
    echo "## Index size:                   $index_size"
else
    echo "## Index size:                   $index_size"
fi

## If a GFF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then
    gff_arg="--sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent"
else
    gff_arg=""
fi

## Overhang length should be read length minus 1
[[ $overhang = "" ]] && overhang=$(( read_len - 1 ))

## Make output dir if needed
mkdir -p "$index_dir"

## Report
echo "## Input FASTA file:             $ref_fa"
[[ "$gff" != "" ]] && echo "## Input GFF file:               $gff"
echo "## Genome index dir (output):    $index_dir"
echo "## Read length:                  $read_len"
echo "## Overhang:                     $overhang"
echo -e "-----------------------\n"


# RUN STAR ---------------------------------------------------------------------
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
echo
