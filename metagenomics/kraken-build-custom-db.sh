#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-kraken-build-%j.out

## Help function
Help() {
    # Display Help
    echo
    echo "## $0: Build a custom Kraken db."
    echo
    echo "## Syntax: $0 -d <kraken-db-dir> -g <genome-dir> [-u <URL>] [-U <URL-file>] [-h]"
    echo 
    echo "## Required options:"
    echo "## -d     Dir for Kraken db (REQUIRED)"
    echo
    echo "## Other options"
    echo "## -g     Dir with/for custom genomes (default: 'genomes' in the Kraken db dir)"
    echo "## -a     Custom genomes are already downloaded - add all genomes in the -g dir"
    echo "## -u     URL for genome FASTA file"
    echo "## -U     File with URLs for genome FASTA files"
    echo "## -h     Print help."
    echo
    echo "## Example: $0 -d refdata/kraken/my-db -g refdata/kraken/my-db/genomes -u https://genome.fa"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}


# SET-UP -----------------------------------------------------------------------
## Report
echo -e "\n## Starting script kraken-build-custom-db.sh..."
date

## Software
### DB-building needs to be done with a manually Kraken2 -- Conda version gives errors
KRAK_BINDIR=/fs/project/PAS0471/jelmer/software/kraken2-2.0.8-beta
KRAK_BIN="$KRAK_BINDIR"/kraken2-build

### Still need to load kraken2-env for "dustmasker", needed to mask low-complexity sequences
### https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#masking-of-low-complexity-sequences
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env

## Bash strict mode
set -euo pipefail

## Option defaults
db_dir=""
genome_dir=""
genome_url=""
genome_url_file=""
genome_dir_all=false

## Get command-line options
while getopts ':d:u:U:g:ah' flag; do
    case "${flag}" in
    d) db_dir="$OPTARG" ;;
    g) genome_dir="$OPTARG" ;;
    u) genome_url="$OPTARG" ;;
    U) genome_url_file="$OPTARG" ;;
    a) genome_dir_all="true" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Derived parameters
lib_dir="$db_dir"/library

## Check input
[[ "$db_dir" = "" ]] && echo -e "\n## ERROR: must specify a db dir with -d" && exit 1
[[ "$genome_url" != "" ]] && [[ "$genome_url_file" != "" ]] && \
    echo -e "\n## ERROR: Use either -u (genome URL) or -U (URL file), not both" && exit 1

## If no URLs were specified, then use all genomes in the genome dir 
[[ "$genome_url" = "" ]] && [[ "$genome_url_file" = "" ]] && \
    genome_dir_all=true && \
    echo -e "\n## NOTE: no URLs specified, using all genomes in genome dir"

## Define genome dir if it wasn't specified with an option
[[ "$genome_dir" = "" ]] && genome_dir="$db_dir"/genomes

## Report
echo
echo "## Database dir:                      $db_dir"
echo "## Genome dir:                        $genome_dir"
echo "## Using all genomes in genome dir:   $genome_dir_all"
[[ "$genome_url" != "" ]] && echo "## Custom genome URL:                 $genome_url"
[[ "$genome_url_file" != "" ]] && echo "## Custom genome URL file:            $genome_url_file"
echo -e "---------------\n\n"

## Make DB directory
mkdir -p "$db_dir"


# DOWNLOAD CUSTOM GENOME(S) ----------------------------------------------------
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


# DOWNLOAD NCBI TAXONOMY -------------------------------------------------------
if [ ! -d "$db_dir"/taxonomy ]; then
    echo -e "\n## Downloading taxonomy..."
    "$KRAK_BIN" --download-taxonomy --db "$db_dir"
fi


# DOWNLOAD TAXON-SPECIFIC LIBRARIES --------------------------------------------
if [ ! -d "$lib_dir"/bacteria ]; then
    echo -e "\n## Downloading library: bacteria..."
    "$KRAK_BIN" --download-library bacteria --db "$db_dir"
fi

if [ ! -d "$lib_dir"/archaea ]; then
    echo -e "\n## Downloading library: archaea..."
    "$KRAK_BIN" --download-library archaea --db "$db_dir"
fi

if [ ! -d "$lib_dir"/viral ]; then
    echo -e "\n## Downloading library: viral..."
    "$KRAK_BIN" --download-library viral --db "$db_dir"
fi

if [ ! -d "$lib_dir"/plasmid ]; then
    echo -e "\n## Downloading library: plasmid..."
    "$KRAK_BIN" --download-library plasmid --db "$db_dir"
fi

if [ ! -d "$lib_dir"/human ]; then
    echo -e "\n## Downloading library: human..."
    "$KRAK_BIN" --download-library human --db "$db_dir"
fi

if [ ! -d "$lib_dir"/UniVec_Core ]; then
    echo -e "\n## Downloading library: UniVec_Core..."
    "$KRAK_BIN" --download-library UniVec_Core --db "$db_dir"
fi


# ADD CUSTOM GENOMES -----------------------------------------------------------
if [ "$genome_url" != "" ]; then
    echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
    "$KRAK_BIN" --add-to-library "$genome_fa" --db "$db_dir"
fi

if [ "$genome_url_file" != "" ]; then
    while read -r genome_url; do
        genome_fa="$genome_dir"/"$(basename "$genome_url")"
        genome_fa=${genome_fa%.*}  # Remove .gz

        echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
        "$KRAK_BIN" --add-to-library "$genome_fa" --db "$db_dir"

    done<"$genome_url_file"
fi

if [ "$genome_dir_all" = true ]; then
    for genome_fa in "$genome_dir"/*; do
        echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
        "$KRAK_BIN" --add-to-library "$genome_fa" --db "$db_dir"
    done
fi


# BUILD THE KRAKEN2 DATABASE ---------------------------------------------------
echo -e "\n## Building final database..."
"$KRAK_BIN" --build \
    --db "$db_dir" \
    --threads "$SLURM_CPUS_ON_NODE"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing file in DB dir:"
ls -lh "$db_dir"
echo -e "\n## Done with script kraken-build-custom-db.sh"
date
