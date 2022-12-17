#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-download-refs-%j.out

## Process command-line arguments
infile=$1        # File with accession numbers for reference genomes, 1 per line
outdir=$2

## Define output file
lookup="$outdir"/lookup.txt

## Software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/entrez-direct
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08

## Report
echo -e "\n## Starting script download-refs.sh"
date
echo "Input file:     $infile"
echo "Output dir:     $outdir"
echo -e "--------------------------\n"

## Make output dir
mkdir -p "$outdir"

echo -e "biosample\tserotype\taccession\tfilename\tfna_presence\tgff_presence" > "$lookup"

## Download genomes
while IFS=$'\t' read -ru9 biosample serotype; do

    echo "$biosample $serotype"
    genomes_file="$outdir"/"$biosample"_"$serotype"_genomeURLs.txt

    esearch -db biosample -query "$biosample" |
        elink -target assembly |
        esummary |
        xtract -pattern DocumentSummary -element FtpPath_GenBank |
        awk -F"/" '{print $0"/"$NF"_genomic.fna.gz" "\n" $0"/"$NF"_genomic.gff.gz" }' > "$genomes_file"
    
    cat "$genomes_file"
    echo
    
    ## Download files
    fna_file="$outdir"/"$(grep "fna.gz" "$genomes_file" | xargs basename | sed 's/\.gz//')"
    gff_file="$outdir"/"$(grep "gff.gz" "$genomes_file" | xargs basename | sed 's/\.gz//')"

    if [[ ! -f "$fna_file" || ! -f "$gff_file" ]]; then
        echo "Downloading files for biosample ${biosample}..."
        wget -q -i "$genomes_file" -P "$outdir"

        ## Unzip
        [[ -f "$fna_file".gz ]] && gunzip -f "$fna_file".gz
        [[ -f "$gff_file".gz ]] && gunzip -f "$gff_file".gz
    else
        echo "Files already present, skipping download for biosample ${biosample}..."
    fi

    ## Check if files were found
    fna_presence=fna_found
    gff_presence=gff_found
    [[ ! -f "$fna_file" ]] && fna_presence=fna_not_found
    [[ ! -f "$gff_file" ]] && gff_presence=gff_not_found

    ## Report to file
    fna_file_base=$(basename "$fna_file")
    accession=$(basename "$fna_file" | sed -E 's/\.[0-9]_.*//')
    echo -e "${biosample}\t${serotype}\t${accession}\t${fna_file_base}\t${fna_presence}\t${gff_presence}" >> "$lookup"

    echo -e "------------------\n"
done 9< "$infile"

## Report
echo -e "\n## Showing contents of lookup file:"
cat "$lookup"
echo -e "\n## Listing files in output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script download-refs.sh"
date
