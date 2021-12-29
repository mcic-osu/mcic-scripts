#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=slurm-tree-build-%j.out

## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2

## Bash strict mode
#set -euo pipefail

## Help
Help() {
    echo
    echo "## $0: Align sequences and build a tree."
    echo
    echo "## Syntax: $0 -i <input-fasta> -o <output-fasta> -t <output-tree> [-p alignment-program] [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     Input FASTA file (REQUIRED)"
    echo "## -o     Output FASTA file (aligned) (REQUIRED)"
    echo "## -t     Output tree file (REQUIRED)"
    echo "## -a     Annotation for tree figure"
    echo "## -f     Addional argument for mafft"
    echo "## -p     Alignment program ('mafft', 'muscle', 'tcoffee', or 'none' (default: 'mafft')"
    echo "## Example: $0 -i in.fa -o out.fa -t out.tree -p muscle"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
align_program="mafft"
annot_in=""
args_mafft=""

## Parse command-line options
while getopts ':i:o:t:p:a:f:h' flag; do
    case "${flag}" in
    i) fa_in="$OPTARG" ;;
    a) annot_in="$OPTARG" ;;
    o) aln="$OPTARG" ;;
    t) tree="$OPTARG" ;;
    p) align_program="$OPTARG" ;;
    f) args_mafft="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process args
outdir_alignment=$(dirname "$aln")
outdir_fasttree=$(dirname "$tree")

tree_fig=${aln%.*}.png

mkdir -p "$outdir_alignment" "$outdir_fasttree"

n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Starting script tree-build.sh..."
date
echo "## Command-line args:"
echo "## Unaligned FASTA file (input):     $fa_in"
echo "## Annotation for tree figure:       $annot_in"
echo "## Aligned FASTA file (output):      $aln"
echo "## Tree file (output):               $tree"
echo "## Tree figure (output):             $tree_fig"
echo "## Alignment program:                $align_program"
[[ "$args_mafft" != "" ]] && echo "## Additional args for mafft:        $args_mafft"
echo -e "--------------------------\n"

## Check input
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input FASTA $fa_in does not exist" && exit 1
[[ "$annot_in" != "" ]] && [[ ! -f "$annot_in" ]] && echo "## ERROR: Input tree annotation file $annot_in does not exist" && exit 1

#? Commands below follow https://www.biorxiv.org/content/10.1101/2020.11.24.396820v1.full
#? GitHub repo: https://github.com/Cyoung02/Phylogenetic-Inference-Benchmarking


# ALIGN ------------------------------------------------------------------------
if [ "$align_program" = mafft ]; then
    echo "## Starting alignment with MAFFT..."
    conda activate mafft-env
    
    mafft --reorder --auto --adjustdirection \
        --leavegappyregion $args_mafft "$fa_in" > "$aln"

elif [ "$align_program" = muscle ]; then
    echo "## Starting alignment with MUSCLE..."
    conda activate muscle-env
    
    muscle -in "$fa_in" -out "$aln"

elif [ "$align_program" = tcoffee ]; then
    echo "## Starting alignment with t-coffee..."
    conda activate t_coffee-env
    
    aln_prefix=${aln%.*}
    
    t_coffee "$fa_in" -n_core="$n_cores" -output=fasta_aln -run_name="$aln_prefix"

    mv "$aln_prefix".fasta_aln "$aln"
else
    echo "## Not running alignment (alignment program variables set to: $align_program)"
fi

## Remove extra info after space from FASTA header lines
## ... and remove "_R_" prefixes for reverse-complemented seqs
sed -i -E -e 's/(^>[^ ]+) .*/\1/' -e 's/^>_R_/>/' "$aln"


# BUILD TREE -------------------------------------------------------------------
echo -e "\n---------------------------"
echo -e "## Starting tree building with FastTree..."
conda activate --stack fasttree-env
cat "$aln" | FastTree -gamma -nt -gtr -out "$tree"


# PLOT TREE --------------------------------------------------------------------
echo -e "\n---------------------------"
echo -e "## Starting tree plotting script..."
module load R/4.0.2-gnu9.1
Rscript scripts/tree-plot.R "$tree" "$aln" "$annot_in" "$tree_fig" #"$show_strain"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n---------------------------"
echo -e "\n## Output files:"
ls -lh "$aln"
ls -lh "$tree"
echo -e "\n## Done with script tree.sh"
date