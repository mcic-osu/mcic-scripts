#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-kraken-build-custom-%j.out

## Help function
Help() {
    echo
    echo "## $0: Build a custom Kraken db."
    echo
    echo "## Syntax: $0 -d <kraken-db-dir> ..."
    echo 
    echo "## Required options:"
    echo "## -d STRING    Dir for Kraken db"
    echo
    echo "## Other options:"
    echo "## -i STRING    Taxonomic groups to include (default: 'abvdhc')"
    echo "                a=archaea, b=bacteria, v=viral, d=plasmid, h=human, c=univec_core, f=fungi, p=plants, z=protozoa"
    echo "                Use a string like 'ab' for archae + bacteria OR 'all' for all groups"
    echo "## -g STRING    Custom genome FASTA file to be added to the db"
    echo "## -G STRING    Dir _with_ (if files are already present) or _for_ (if files are to be downloaded) custom genome FASTA files"
    echo "## -u STRING    URL for genome FASTA file to be downloaded and added to the db"
    echo "## -U STRING    File with URLs for genome FASTA files to be downloaded and added to the db"
    echo "## -h           Print this help message"
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
### DB-building needs to be done with a manually installed Kraken2 -- Conda version gives errors
#KRAK_BINDIR=/fs/project/PAS0471/jelmer/software/kraken2-2.0.8-beta
#KRAK_BIN="$KRAK_BINDIR"/kraken2-build

### Still need to load kraken2-env for "dustmasker", needed to mask low-complexity sequences
### https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#masking-of-low-complexity-sequences
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env

## Bash strict mode
set -euo pipefail

## Option defaults
db_dir=""
genome_dir=""
genome_url=""
genome_url_file=""
genome_dir_all=false
genome_fa=""
include_taxa="abvdhc"

include_archaea=false
include_viral=false
include_bacteria=false
include_plasmid=false
include_human=false
include_univec=false
include_fungi=false
include_plants=false
include_protozoa=false

## Get command-line options
while getopts ':d:u:U:g:G:i:h' flag; do
    case "${flag}" in
    d) db_dir="$OPTARG" ;;
    g) genome_fa="$OPTARG" ;;
    G) genome_dir="$OPTARG" ;;
    u) genome_url="$OPTARG" ;;
    U) genome_url_file="$OPTARG" ;;
    i) include_taxa="$OPTARG" ;;
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

## If no URLs/FASTA files were specified, but a genome dir was specified, then use all genomes in the genome dir 
[[ "$genome_url" = "" ]] && [[ "$genome_url_file" = "" ]] && [[ "$genome_fa" = "" ]] && [[ "$genome_dir" != "" ]] && \
    genome_dir_all=true && \
    echo -e "\n## NOTE: Using all genomes in genome dir $genome_dir"

## Figure out which taxa to use
if [[ "$include_taxa" = "all" ]]; then
    include_archaea=true
    include_viral=true
    include_bacteria=true
    include_plasmid=true
    include_human=true
    include_univec=true
    include_fungi=true
    include_plants=true
    include_protozoa=true
else
    [[ "$include_taxa" = *a* ]] && include_archaea=true
    [[ "$include_taxa" = *b* ]] && include_viral=true
    [[ "$include_taxa" = *v* ]] && include_bacteria=true
    [[ "$include_taxa" = *d* ]] && include_plasmid=true
    [[ "$include_taxa" = *h* ]] && include_human=true
    [[ "$include_taxa" = *c* ]] && include_univec=true
    [[ "$include_taxa" = *f* ]] && include_fungi=true
    [[ "$include_taxa" = *p* ]] && include_plants=true
    [[ "$include_taxa" = *z* ]] && include_protozoa=true
fi

## Report
echo
echo "## Database dir:                 $db_dir"
[[ "$genome_dir" != "" ]] && echo "## Genome dir:                   $genome_dir"
[[ "$genome_fa" != "" ]] && echo "## Custom genome FASTA:          $genome_fa"
[[ "$genome_url" != "" ]] && echo "## Custom genome URL:            $genome_url"
[[ "$genome_url_file" != "" ]] && echo "## Custom genome URL file:       $genome_url_file"
echo
echo "## Including taxa:"
[[ "$include_archaea" = true ]] && echo "archaea"
[[ "$include_viral" = true ]] && echo "viral"
[[ "$include_bacteria" = true ]] && echo "bacteria"
[[ "$include_plasmid" = true ]] && echo "plasmid"
[[ "$include_human" = true ]] && echo "human"
[[ "$include_univec" = true ]] && echo "univec"
[[ "$include_fungi" = true ]] && echo "fungi"
[[ "$include_plants" = true ]] && echo "plants"
[[ "$include_protozoa" = true ]] && echo "protozoa"
echo -e "---------------\n"

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
        genome_fa_one="$genome_dir"/"$(basename "$genome_url")"
        genome_fa_one=${genome_fa_one%.*}  # Remove .gz

        echo "## Custom genome: $genome_fa_one"

        if [ ! -f "$genome_fa_one" ]; then
            echo "## Custom genome not present, will download..."
            wget -q -P "$genome_dir" "$genome_url"
            echo "## Unzipping custom genome..."
            gunzip "$genome_fa_one".gz
        else
            echo "## Custom genome already present."
        fi

        echo "## Listing custom genome:"
        ls -lh "$genome_fa_one"
        echo
    done<"$genome_url_file"
fi


# DOWNLOAD NCBI TAXONOMY -------------------------------------------------------
if [ ! -d "$db_dir"/taxonomy ]; then
    echo -e "\n## Downloading taxonomy..."
    kraken2-build --download-taxonomy --db "$db_dir"
fi


# DOWNLOAD TAXON-SPECIFIC LIBRARIES --------------------------------------------
if [ ! -d "$lib_dir"/archaea ] && [ "$include_archaea" = true ]; then
    echo -e "\n## Downloading library: archaea..."
    kraken2-build --download-library archaea --db "$db_dir"
fi

if [ ! -d "$lib_dir"/bacteria ] && [ "$include_bacteria" = true ]; then
    echo -e "\n## Downloading library: bacteria..."
    kraken2-build --download-library bacteria --db "$db_dir"
fi

if [ ! -d "$lib_dir"/plasmid ] && [ "$include_plasmid" = true ]; then
    echo -e "\n## Downloading library: plasmid..."
    kraken2-build --download-library plasmid --db "$db_dir"
fi

if [ ! -d "$lib_dir"/viral ] && [ "$include_viral" = true ]; then
    echo -e "\n## Downloading library: viral..."
    kraken2-build --download-library viral --db "$db_dir"
fi

if [ ! -d "$lib_dir"/human ] && [ "$include_human" = true ]; then
    echo -e "\n## Downloading library: human..."
    kraken2-build --download-library human --db "$db_dir"
fi

if [ ! -d "$lib_dir"/fungi ] && [ "$include_fungi" = true ]; then
    echo -e "\n## Downloading library: fungi..."
    kraken2-build --download-library fungi --db "$db_dir"
fi

if [ ! -d "$lib_dir"/plant ] && [ "$include_archaea" = true ]; then
    echo -e "\n## Downloading library: plant..."
    kraken2-build --download-library plant --db "$db_dir"
fi

if [ ! -d "$lib_dir"/protozoa ] && [ "$include_protozoa" = true ]; then
    echo -e "\n## Downloading library: protozoa..."
    kraken2-build --download-library protozoa --db "$db_dir"
fi

if [ ! -d "$lib_dir"/UniVec_Core ] && [ "$include_univec" = true ]; then
    echo -e "\n## Downloading library: UniVec_Core..."
    kraken2-build --download-library UniVec_Core --db "$db_dir"
fi


# ADD CUSTOM GENOMES -----------------------------------------------------------
if [ "$genome_fa" != "" ]; then
    
    ## If there is a single genome to be added
    echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
    kraken2-build --add-to-library "$genome_fa" --db "$db_dir"

elif [ "$genome_dir_all" = true ]; then
    
    ## If all genomes in a dir should be added
    for genome_fa in "$genome_dir"/*; do
        echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
        kraken2-build --add-to-library "$genome_fa" --db "$db_dir"
    done

elif [ "$genome_url_file" != "" ]; then
    
    ## If a file with URLs was provided
    while read -r genome_url; do
        genome_fa="$genome_dir"/"$(basename "$genome_url")"
        genome_fa=${genome_fa%.*}  # Remove .gz

        echo -e "\n## Adding custom genome $genome_fa to Kraken library..."
        kraken2-build --add-to-library "$genome_fa" --db "$db_dir"
    done<"$genome_url_file"

else
    echo -e "\n## Not adding any custom genomes..."
fi


# BUILD THE KRAKEN2 DATABASE ---------------------------------------------------
echo -e "\n## Building the final database..."
kraken2-build \
    --build \
    --db "$db_dir" \
    --threads "$SLURM_CPUS_ON_NODE"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing file in DB dir:"
ls -lh "$db_dir"
echo -e "\n## Done with script kraken-build-custom-db.sh"
date
