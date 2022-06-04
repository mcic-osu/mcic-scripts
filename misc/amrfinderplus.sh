#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=amrfinderplus
#SBATCH --output=slurm-amrfinderplus-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run AMRFinderPlus to identify antimicrobial resistance genes and SNPs in a microbial genome assembly"
  echo
  echo "Syntax: $0 -n <input-fna> -a <input-faa> -g <input-gff> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -n STRING         Nucleotide FASTA input file (genome assembly)"
  echo "    -a STRING         Amino acid (protein) FASTA input file (proteome)"
  echo "    -g STRING         GFF input file (annotation)"
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -s                Organism                      [default: none]"
  echo "                          For a list of options, run 'amrfinder -l'"
  echo "                          See https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option"
  echo "    -x                Don't change GFF file"
  echo "                          The default is to change 'ID=' entries to 'Name='."
  echo "                          This is necessary at least if the GFF was produced by Prokka"
  echo "                          See https://github.com/ncbi/amr/issues/26"
  echo "    -a STRING         Other argument(s) to pass to AMRFinderPlus"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:    $0 -n results/spades/smpA.fasta -a results/prokka/smpA.faa -g results/prokka/smpA.gff  -o results/ksnp3"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "AMRFinderPlus documentation: https://github.com/ncbi/amr/wiki"
  echo "AMRFinderPlus paper: https://www.nature.com/articles/s41598-021-91456-0"
  echo
}

## Option defaults
assembly_fna=""
assembly_faa=""
gff=""
outdir=""
organism=""
change_gff="true"
more_args=""

## Parse command-line options
while getopts ':n:p:g:o:s:a:hx' flag; do
  case "${flag}" in
    n) assembly_fna="$OPTARG" ;;
    p) assembly_faa="$OPTARG" ;;
    g) gff="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    s) organism="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    x) change_gff="false" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Check input
[[ "$assembly_fna" = "" ]] && echo "## ERROR: Please specify a nucleotide FASTA input file with -n" >&2 && exit 1
[[ "$assembly_faa" = "" ]] && echo "## ERROR: Please specify a protein FASTA input file with -p" >&2 && exit 1
[[ "$gff" = "" ]] && echo "## ERROR: Please specify a GFF input file with -g" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1

[[ ! -f "$assembly_fna" ]] && echo "## ERROR: Input file (-n) $assembly_fna does not exist" >&2 && exit 1
[[ ! -f "$assembly_faa" ]] && echo "## ERROR: Input file (-p) $assembly_faa does not exist" >&2 && exit 1
[[ ! -f "$gff" ]] && echo "## ERROR: Input file (-g) $gff does not exist" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/amrfinderplus-3.10.30

## Bash strict settings
set -euo pipefail

## Define output files etc
sampleID=$(basename "$assembly_fna" .fasta)
outfile="$outdir"/"$sampleID".txt
mutation_report="$outdir"/"$sampleID"_mutation-report.txt

## Build organism argument
if [[ $organism != "" ]]; then
    organism_arg="--organism $organism"
else
    organism_arg=""
fi

## Report
echo "## Starting script amrfinderplus.sh"
date
echo
echo "## Input nucleotide FASTA file:          $assembly_fna"
echo "## Input protein FASTA file:             $assembly_faa"
echo "## Input GFF file:                       $gff"
echo "## Output dir:                           $outdir"
[[ $organism != "" ]] && echo "## Organism:                             $organism"
[[ $more_args != "" ]] && echo "## Other arguments to pass to AmrFinderPlus:    $more_args"
echo -e "--------------------\n"

## Change GFF to be compliant with Amrfinderplus
if [[ "$change_gff" = true ]]; then
    echo "## Now editing the GFF file..."    
    sed -E 's/Name=[^;]+;//' "$gff" | sed 's/ID=/Name=/' > "$outdir"/"$sampleID".gff
    gff="$outdir"/"$sampleID".gff
fi

## Make output dir
mkdir -p "$outdir"


# RUN AmrfinderPlus --------------------------------------------------------------------
echo -e "\n## Now running AmrfinderPlus...."
amrfinder \
    --nucleotide "$assembly_fna" \
    --protein "$assembly_faa" \
    --gff "$gff" \
    --threads "$SLURM_CPUS_PER_TASK" \
    -o "$outfile" \
    --mutation_all "$mutation_report" \
    --name "$sampleID" $organism_arg $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"/"$sampleID"*
echo -e "\n## Done with script amrfinderplus.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
