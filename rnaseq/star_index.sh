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
fasta=""
outdir=""
index_size="auto"
gff=""
overhang="auto"
read_len=150

## Parse command-line options
while getopts ':i:a:o:s:r:v:h' flag; do
    case "${flag}" in
        i) fasta="$OPTARG" ;;
        a) gff="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
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
source activate /fs/project/PAS0471/jelmer/conda/star-2.7.10a

## Strict bash settings
set -euo pipefail
## Check input
[[ "$fasta" = "" ]] && Die "Please specify an input file with -i/--infile" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$fasta" ]] && Die "Input FASTA file $fasta does not exist"
[[ "$gff" != "" ]] && [[ ! -f "$gff" ]] && Die "Input GFF file $fasta does not exist"


## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT STAR_INDEX.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $fasta"
[[ "$gff" != "" ]] && echo "## Input GFF file:                   $gff"
echo "Output dir:                       $outdir"
echo "Read length:                      $read_len"
echo "Memory in GB / bytes:             $mem_gb / $mem_bytes"
[[ $more_args != "" ]] && echo "Other arguments for STAR:         $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$fasta"
echo
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

## STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $fasta = *gz ]]; then
    fasta_unzip=${fasta/.gz/}
    if [[ ! -f $fasta_unzip ]]; then
        echo "Unzipping gzipped FASTA file..."
        gunzip -c "$fasta" > "$fasta_unzip"
    else
        echo "Unzipped version of the FASTA file already exists: $fasta_unzip"
    fi
    fasta="$fasta_unzip"
fi

## Determine index size
if [ "$index_size" = "auto" ]; then
    echo -e "\nAutomatically determining the index size..."
    genome_size=$(grep -v "^>" "$fasta" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    echo "Genome size (autom. determined):  $genome_size"
    echo "Index size (autom. determined):   $index_size"
else
    echo "Index size:                       $index_size"
fi

## If a GFF file is provided, build the appropriate argument for STAR
if [ "$gff" != "" ]; then
    ## Overhang length should be read length minus 1 - only if GFF is included
    [[ $overhang = "auto" ]] && overhang=$(( read_len - 1 ))
    gff_arg="--sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $overhang"
    echo "Overhang:                     $overhang"
else
    gff_arg=""
fi

## Report
echo "=========================================================================="
## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


## Run STAR
echo -e "\n# Indexing with STAR...."
${e}Time STAR \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM "$mem_bytes" \
    --genomeDir "$outdir" \
    --genomeFastaFiles "$fasta" \
    --genomeSAindexNbases "$index_size" \
    --runThreadN "$threads" ${gff_arg}


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
