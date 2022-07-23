#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=STAR_index
#SBATCH --output=slurm-STAR-index-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Index a reference genome FASTA file with STAR."
  echo
  echo "Syntax: $0 -i <input-FASTA> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "   -i FILE        Input reference FASTA file"
  echo "   -o DIR         Output directory for index files (will be created if needed)"
  echo
  echo "Other options:"
  echo "   -a FILE        Reference annotation (GFF/GTF) file [default: no GFF/GTF, but this is not recommended]"
  echo "   -s INTEGER     Index size                          [default: 'auto' => automatically determined from genome size]"
  echo "   -r INTEGER     Read length                         [default: '150' (bp)]"
  echo "   -v INTEGER     Overhang                            [default: 'auto' => read length minus 1]"
  echo "                  (Note: overhang only applies if GFF/GTF file is provided!)"
  echo "   -h             Print this help message and exit"
  echo
  echo "Example:       $0 -i refdata/my_genome.fa -o refdata/star_index -a refdata/my_genome.gff"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Note: the script will check how much memory has been allocated to the SLURM job (default: 64GB),"
  echo "      and pass that to STAR via the 'limitGenomeGenerateRAM argument'."
  echo "      When allocating more memory to the SLURM job,"
  echo "      wich can be necessary for large genomes, this will therefore be passed to STAR as well."
  echo
}

## Option defaults
ref_fa=""
index_dir=""
index_size="auto"
gff=""
overhang="auto"
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
## Check inputs
[[ ! -f "$ref_fa" ]] && echo "## ERROR: Input FASTA (-i) $ref_fa does not exist" >&2 && exit 1
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && echo "## ERROR: Input file GFF (-a) $gff does not exist" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/.conda/envs/star-env

## Strict bash settings
set -euo pipefail

## Determine amount of memory
mem_bytes=$((SLURM_MEM_PER_NODE * 1000000))
mem_gb=$((SLURM_MEM_PER_NODE / 1000))


## Report
echo
echo "## Starting script star_index.sh"
date
echo
echo "## Input FASTA file:                 $ref_fa"
[[ "$gff" != "" ]] && echo "## Input GFF file:                   $gff"
echo "## Genome index dir (output):        $index_dir"
echo "## Read length:                      $read_len"
echo "## Memory in GB / bytes:             $mem_gb / $mem_bytes"
echo

## Make output dir if needed
mkdir -p "$index_dir"


# PREP INPUT -------------------------------------------------------------------
## STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $ref_fa = *gz ]]; then
    ref_fa_unzip=${ref_fa/.gz/}
    if [[ ! -f $ref_fa_unzip ]]; then
        echo "## Unzipping gzipped FASTA file..."
        gunzip -c "$ref_fa" > "$ref_fa_unzip"
    else
        echo "## Unzipped version of the FASTA file already exists: $ref_fa_unzip"
    fi
    ref_fa="$ref_fa_unzip"
fi

## Determine index size
if [ "$index_size" = "auto" ]; then
    echo -e "\n## Automatically determining the index size..."
    genome_size=$(grep -v "^>" "$ref_fa" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    echo "## Genome size (autom. determined):  $genome_size"
    echo "## Index size (autom. determined):   $index_size"
else
    echo "## Index size:                       $index_size"
fi

## If a GFF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then
    ## Overhang length should be read length minus 1 - only if GFF is included
    [[ $overhang = "auto" ]] && overhang=$(( read_len - 1 ))
    gff_arg="--sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $overhang"
    echo "## Overhang:                     $overhang"
else
    gff_arg=""
fi

echo -e "-----------------------\n"


# RUN STAR ---------------------------------------------------------------------
echo -e "\n## Indexing genome with STAR...."
STAR --runMode genomeGenerate \
     --limitGenomeGenerateRAM "$mem_bytes" \
     --genomeDir "$index_dir" \
     --genomeFastaFiles "$ref_fa" \
     --genomeSAindexNbases "$index_size" \
     --runThreadN "$SLURM_CPUS_ON_NODE" ${gff_arg}


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$index_dir"
echo -e "\n## Done with script star_index.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
