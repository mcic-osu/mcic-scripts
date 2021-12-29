#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-kraken-build-%j.out

## Software
### DB-building needs to be done with a manually Kraken2 -- Conda version gives errors
KRAKEN_BIN_DIR=software/kraken/kraken2-2.0.8-beta

### Still need to load kraken2-env for "dustmasker", needed to mask low-complexity sequences
### https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#masking-of-low-complexity-sequences
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate kraken2-env

## Bash strict mode
set -euo pipefail

## Help
Help() {
    # Display Help
    echo
    echo "## $0: Build a custom Kraken db."
    echo
    echo "## Syntax: $0 -d <kraken-db-dir> -g <genome-dir> -f <forward-primer> -r <reverse-primer> [-d] [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -d     Kraken db dir (REQUIRED)"
    echo "## -g     Genome dir (REQUIRED)"
    echo "## -u     URL for genome FASTA file"
    echo "## -U     File with URLs for genome FASTA files"
    echo "## Example: $0 -i refdata/kraken/my-db -u https://genome.fa -g refdata/kraken/my-db/genomes"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
db_dir=""
genome_dir=""
genome_url=""
genome_url_file=""
genome_dir_all=false

## Get command-line options
while getopts ':d:u:U:g:h' flag; do
    case "${flag}" in
    d) db_dir="$OPTARG" ;;
    g) genome_dir="$OPTARG" ;;
    u) genome_url="$OPTARG" ;;
    U) genome_url_file="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process options
[[ "$db_dir" = "" ]] && echo "ERROR: must specift a db dir with -d" && exit 1
[[ "$genome_url" != "" ]] && [[ "$genome_url_file" != "" ]] && echo "ERROR: Use either -u (genome URL) or -U (URL file), not both" && exit 1
[[ "$genome_url" = "" ]] && [[ "$genome_url_file" = "" ]] && genome_dir_all=true && echo "NOTE: no URLs specified, using all genomes in genome dir" && exit 1
[[ "$genome_dir" = "" ]] && genome_dir="$db_dir"/genomes

## Report
echo "## Starting script kraken-build-custom-db.sh..."
date
echo "## Database dir: $db_dir"
echo "## Genome dir: $genome_dir"
[[ "$genome_url" != "" ]] && echo "## Custom genome URL: $genome_url"
[[ "$genome_url_file" != "" ]] && echo "## Custom genome URL file: $genome_url_file"
echo "## Using all genomes in genome dir: $genome_dir_all"
echo -e "---------------\n\n"

## Make DB directory
mkdir -p "$db_dir"

## Download custom genome(s)
if [ "$genome_url" != "" ]; then
    genome_fa="$genome_dir"/"$(basename "$genome_url")"
    genome_fa=${genome_fa%.*}  # Remove .gz

    if [ ! -f "$genome_fa" ]; then
        echo "## Custom genome not present, will download..."
        wget -q -P "$genome_dir" "$genome_url"
        echo "## Unzipping custom genome..."
        gunzip "$genome_fa".gz
    else
        echo "## Custom genome already present."
    fi

    echo "## Listing custom genome:"
    ls -lh "$genome_fa"
fi

if [ "$genome_url_file" != "" ]; then
    
    while read -r genome_url; do
        genome_fa="$genome_dir"/"$(basename "$genome_url")"
        genome_fa=${genome_fa%.*}  # Remove .gz

        echo "## Custom genome: $genome_fa"

        if [ ! -f "$genome_fa" ]; then
            echo "## Custom genome not present, will download..."
            wget -q -P "$genome_dir" "$genome_url"
            echo "## Unzipping custom genome..."
            gunzip "$genome_fa".gz
        else
            echo "## Custom genome already present."
        fi

        echo "## Listing custom genome:"
        ls -lh "$genome_fa"
        echo
    
    done<"$genome_url_file"

fi

## Download taxonomy
if [ ! -d "$db_dir"/taxonomy ]; then
    echo -e "\n## Downloading taxonomy..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-taxonomy --db "$db_dir"
fi

## Download libraries
lib_dir="$db_dir"/library

if [ ! -d "$lib_dir"/bacteria ]; then
    echo -e "\n## Downloading library: bacteria..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library bacteria --db "$db_dir"
fi

if [ ! -d "$lib_dir"/archaea ]; then
    echo -e "\n## Downloading library: archaea..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library archaea --db "$db_dir"
fi

if [ ! -d "$lib_dir"/viral ]; then
    echo -e "\n## Downloading library: viral..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library viral --db "$db_dir"
fi

if [ ! -d "$lib_dir"/plasmid ]; then
    echo -e "\n## Downloading library: plasmid..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library plasmid --db "$db_dir"
fi

if [ ! -d "$lib_dir"/human ]; then
    echo -e "\n## Downloading library: human..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library human --db "$db_dir"
fi

if [ ! -d "$lib_dir"/UniVec_Core ]; then
    echo -e "\n## Downloading library: UniVec_Core..."
    "$KRAKEN_BIN_DIR"/kraken2-build --download-library UniVec_Core --db "$db_dir"
fi

## Add custom genome
if [ "$genome_url" != "" ]; then
    echo -e "\n## Adding custom genome..."
    "$KRAKEN_BIN_DIR"/kraken2-build --add-to-library "$genome_fa" --db "$db_dir"
fi

if [ "$genome_url_file" != "" ]; then
    
    while read -r genome_url; do
        genome_fa="$genome_dir"/"$(basename "$genome_url")"
        genome_fa=${genome_fa%.*}  # Remove .gz

        echo -e "\n## Adding custom genome $genome_fa..."
        "$KRAKEN_BIN_DIR"/kraken2-build --add-to-library "$genome_fa" --db "$db_dir"

    done<"$genome_url_file"

fi

if [ "$genome_dir_all" = true ]; then
    
    for genome_fa in "$genome_dir"/*; do
        echo -e "\n## Adding custom genome $genome_fa..."
        "$KRAKEN_BIN_DIR"/kraken2-build --add-to-library "$genome_fa" --db "$db_dir"
    done

fi

## Build database
echo -e "\n## Building final database..."
"$KRAKEN_BIN_DIR"/kraken2-build --build --threads "$SLURM_CPUS_ON_NODE" --db "$db_dir"

## Report:
echo -e "\n## Done with script kraken-build-custom-db.sh"
date