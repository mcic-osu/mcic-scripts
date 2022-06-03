#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=amrfinderplus
#SBATCH --output=slurm-amrfinderplus-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run amrfinderplus."
  echo
  echo "Syntax: $0 -i <input-file> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "-i STRING         Input file which should contain one row per genome and two columns:"
  echo "                     - The first column has the path to a genome FASTA"
  echo "                     - The second column has an ID for that genome"
  echo "                  E.g., a file for two genomes could look like this:"
  echo "                  results/spades/sampleA/contigs.fasta sampleA"
  echo "                  results/spades/sampleB/contigs.fasta sampleB"
  echo "-o STRING         Output dir"
  echo
  echo "Other options:"
  echo "-k INTEGER        K-mer size (odd integer)          [default: automatically determined]"
  echo "                  If you don't provide a k-mer size, the script will automatically"
  echo "                     determine one using the kSNP3 utility program Kchooser"
  echo "-a STRING         Other argument(s) to pass to kSNP3"
  echo "                  Note that by default, the script will run kSNP3 with the following optional arguments:"
  echo "                     -vcf      To also output a VCF file"
  echo "                     -ML       To also create an ML tree"
  echo "                     -core     To also output files for 'core SNPs' only"
  echo "-h                Print this help message and exit"
  echo
  echo "Example: $0 -i infile_list.txt -o results/ksnp3"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "kSNP3 documentation (PDF download): https://sourceforge.net/projects/ksnp/files/kSNP3.1.2%20User%20Guide%20.pdf/download"
  echo "kSNP3 paper: https://academic.oup.com/bioinformatics/article/31/17/2877/183216"
  echo
}

## Option defaults
infile=""
outdir=""
organism="Salmonella"
more_args=""

## Parse command-line options
while getopts ':i:o:O:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    O) organism="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Check input
[[ "$infile" = "" ]] && echo "## ERROR: Please specify an input dir with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file (-i) $infile does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/amrfinderplus-3.10.30

## Bash strict settings
set -euo pipefail

## Define output files etc
sampleID=$(basename "$infile" .fasta)
outfile="$outdir"/"$sampleID".txt
mutation_report="$outdir"/"$sampleID"_mutation-report.txt

## Report
echo "## Starting script amrfinderplus.sh"
date
echo
echo "## Input file:                           $infile"
echo "## Output dir:                           $outdir"
echo "## Organism:                             $organism"
[[ $more_args != "" ]] && echo "## Other arguments to pass to AmrFinderPlus:    $more_args"
echo -e "--------------------\n"

## Make output dir
mkdir -p "$outdir"


# RUN AmrfinderPlus --------------------------------------------------------------------
echo -e "\n## Now running AmrfinderPlus...."
amrfinder \
    --nucleotide "$infile" \
    --organism "$organism" \
    --threads "$SLURM_CPUS_PER_TASK" \
    -o "$outfile" \
    --mutation_all "$mutation_report" \
    --name "$sampleID" $more_args

# --protein <protein_fasta>
# --gff <gff_file>


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script amrfinderplus.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
