#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=slurm-tree-build-%j.out

## Help
Help() {
    echo
    echo "## $0: Align sequences and build a tree."
    echo
    echo "## Syntax: $0 -i <input-fasta> -o <output-fasta> -t <output-tree> [-p alignment-program] [-h]"
    echo "## Options:"
    echo "## -i STR     Input FASTA file (REQUIRED)"
    echo "## -o STR     Output (aligned) FASTA filename (REQUIRED)"
    echo "## -t STR     Output tree filename (REQUIRED)"
    echo
    echo "## -p STR     Alignment program ('mafft', 'muscle', 'tcoffee', or 'none' (default: 'mafft')"
    echo "## -f STR     Additional argument(s) for mafft"
    echo "## -l         Make a plot of the tree (default: no plot)"
    echo "## -a STR     Annotation df for tree plot"
    echo
    echo "## -h         Print help."
    echo
    echo "## Example: $0 -i in.fa -o out.fa -t out.tree -p muscle"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Software and scripts
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
conda activate fasttree-env
#? Note: additional Conda environments are loaded below, depending on the alignment program, etc

## Bash strict mode
set -euo pipefail

## Option defaults
align_program="mafft"
annot_in=""
args_mafft=""

## Parse command-line options
while getopts ':i:o:t:p:a:f:lh' flag; do
    case "${flag}" in
    i) fa_in="$OPTARG" ;;
    a) annot_in="$OPTARG" ;;
    o) aln="$OPTARG" ;;
    t) tree="$OPTARG" ;;
    p) align_program="$OPTARG" ;;
    f) args_mafft="$OPTARG" ;;
    l) make_plot=true ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process parameters
### Output dirs
outdir_alignment=$(dirname "$aln")
outdir_fasttree=$(dirname "$tree")
mkdir -p "$outdir_alignment" "$outdir_fasttree"
### Tree figure
tree_fig=${aln%.*}.png
### Number of cores
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Starting script tree-build.sh..."
date
echo "## Unaligned FASTA file (input):     $fa_in"
echo "## Aligned FASTA file (output):      $aln"
echo "## Tree file (output):               $tree"
echo "## Alignment program:                $align_program"
[[ "$make_plot" = true ]] && echo "## Tree figure (output):             $tree_fig"
[[ "$annot_in" != "" ]] && echo "## Annotation for tree figure:       $annot_in"
[[ "$args_mafft" != "" ]] && echo "## Additional args for mafft:        $args_mafft"
echo -e "--------------------------\n"

## Check input
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input FASTA $fa_in does not exist" && exit 1
[[ "$annot_in" != "" ]] && [[ ! -f "$annot_in" ]] && echo "## ERROR: Input tree annotation file $annot_in does not exist" && exit 1


# ALIGN ------------------------------------------------------------------------
if [ "$align_program" = mafft ]; then
    echo "## Starting alignment with MAFFT..."
    conda activate --stack mafft-env
    
    mafft --reorder --auto --adjustdirection \
        --leavegappyregion $args_mafft "$fa_in" > "$aln"

elif [ "$align_program" = muscle ]; then
    echo "## Starting alignment with MUSCLE..."
    conda activate --stack muscle-env
    
    muscle -in "$fa_in" -out "$aln"

elif [ "$align_program" = tcoffee ]; then
    echo "## Starting alignment with t-coffee..."
    conda activate --stack t_coffee-env
    
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
cat "$aln" | FastTree -gamma -nt -gtr -out "$tree"


# PLOT TREE --------------------------------------------------------------------
if [ "$make_plot" = true ]; then
    echo -e "\n---------------------------"
    echo -e "## Starting tree plotting script..."

    module load R/4.0.2-gnu9.1
    TREE_PLOT_SCRIPT=mcic-scripts/trees/tree-plot.R

    if [ "$annot_in" != "" ]; then
        Rscript "$TREE_PLOT_SCRIPT" -t "$tree" -a "$aln" -o "$tree_fig" -n "$annot_in"
    else
        Rscript "$TREE_PLOT_SCRIPT" -t "$tree" -a "$aln" -o "$tree_fig"
    fi
fi

# WRAP UP ----------------------------------------------------------------------
echo -e "\n---------------------------"
echo -e "\n## Output files:"
ls -lh "$aln"
ls -lh "$tree"
ls -lh "$tree_fig"
echo -e "\n## Done with script tree.sh"
date


#? Commands follow https://www.biorxiv.org/content/10.1101/2020.11.24.396820v1.full
#? GitHub repo: https://github.com/Cyoung02/Phylogenetic-Inference-Benchmarking