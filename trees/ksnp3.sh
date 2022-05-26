#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=ksnp3
#SBATCH --output=slurm-ksnp3-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run kSNP3 to identity SNPs among a set of bacterial genomes."
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
kmer_size="auto"
more_args=""

## Parse command-line options
while getopts ':i:o:k:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    k) kmer_size="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## If needed, make paths absolute because we have to move into the outdir
[[ ! $outdir =~ ^/ ]] && outdir="$PWD"/"$outdir"
[[ ! $infile =~ ^/ ]] && outdir="$PWD"/"$infile"

## Check input
[[ "$infile" = "" ]] && echo "## ERROR: Please specify an input dir with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file (-i) $infile does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/knsp-3.1

## Bash strict settings
set -euo pipefail

## Report
echo "## Starting script ksnp3.sh"
date
echo
echo "## Input file:                           $infile"
echo "## Output dir:                           $outdir"
echo "## Kmer size:                            $kmer_size"
[[ $more_args != "" ]] && echo "## Other arguments to pass to kNSP3:    $more_args"
echo "## Showing contents of input file with file names:"
cat "$infile"
echo -e "--------------------\n"

## Make output dir
mkdir -p "$outdir"


# RUN kSNP3 --------------------------------------------------------------------
## Move into the output dir because kSNP3 will dump files into the working dir
cd "$outdir" || exit 1

## Determine kmer-size
if [ "$kmer_size" = "auto" ]; then
    echo -e "\n## Because kmer size is set to 'auto', will now pick a kmer size..."
    
    ## Create combined FASTA file
    echo "## Creating a merged FASTA file for Kchooser..."
    MakeFasta "$infile" "$outdir"/merged_fa_for_kchooser.fa
    
    ## Run Kchooser
    echo -e "\n## Running Kchooser to pick a kmer size..."
    Kchooser "$outdir"/merged_fa_for_kchooser.fa
    
    ## Get the optimum kmer size from the Kchooser report
    if grep -q "The optimum value of K is" Kchooser.report; then
        kmer_size=$(grep "The optimum value of K is" Kchooser.report | sed -E 's/.*is ([0-9]+)./\1/')
    else
        echo -e "\n## ERROR: No optimum K-value found in Kchooser.report" >&2 && exit 1
    fi
    
    ## Report
    echo -e "\n## The optimum kmer size is:         $kmer_size"

    ## A kmer size of 31 probably means Kchooser needs to be rerun
    if [ "$kmer_size" = 31 ]; then
        echo
        echo "## ERROR: The optimal kmer-size is 31, which is probably incorrect." >&2
        echo "##        See the kSNP3 documentation to rerun Kchooser with a different fraction of unique kmers" >&2 && exit 1
    fi
fi

## Run kSNP3
echo -e "\n## Now running kSNP3...."
kSNP3 \
    -in "$infile" \
    -outdir "$outdir" \
    -k "$kmer_size" \
    -vcf -ML \
    -CPU "$SLURM_CPUS_PER_TASK"

#? -ML     Also output an ML tree
#? -vcf    Also output a VCF file


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script ksnp3.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
