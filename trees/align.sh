#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --output=slurm-align-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Align sequences in a multi-FASTA file"
    echo
    echo "Syntax: $0 -i <input-fasta> -o <output-fasta>"
    echo
    echo "Required options:"
    echo "  -i STRING       Input (unaligned) FASTA file"
    echo "  -o STR          Output (aligned) FASTA file"
    echo
    echo "Other options:"
    echo "  -p STRING       Alignment program ('mafft', 'muscle', 'tcoffee', or 'none') [default: 'mafft']"
    echo "  -f STRING       Additional argument(s) for mafft"
    echo "  -h              Print this help message and exit"
    echo
    echo "## Example:       $0 -i in.fa -o aligned.fa -p muscle"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
fa_in=""
fa_out=""
align_program="mafft"
args_mafft=""

## Parse command-line options
while getopts ':i:o:p:f:h' flag; do
    case "${flag}" in
    i) fa_in="$OPTARG" ;;
    o) fa_out="$OPTARG" ;;
    p) align_program="$OPTARG" ;;
    f) args_mafft="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Bash strict mode
set -euo pipefail

## Software and scripts
module load python/3.6-conda5.2
#? Note: Conda environments are loaded below, depending on the alignment program, etc

## Check input
[[ "$fa_in" = "" ]] && echo "## ERROR: Please specify a FASTA input file with -i " >&2 && exit 1
[[ "$fa_out" = "" ]] && echo "## ERROR: Please specify a FASTA output file with -o " >&2 && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input FASTA $fa_in does not exist" >&2 && exit 1

## Process parameters - output dirs
outdir_alignment=$(dirname "$fa_out")
mkdir -p "$outdir_alignment"

## Process parameters - number of cores
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo
echo "## Starting script tree-build.sh..."
date
echo
echo "## Unaligned FASTA file (input):     $fa_in"
echo "## Aligned FASTA file (output):      $fa_out"
echo "## Alignment program:                $align_program"
[[ "$args_mafft" != "" ]] && echo "## Additional args for mafft:        $args_mafft"
echo -e "--------------------------\n"


# ALIGN ------------------------------------------------------------------------
if [ "$align_program" = mafft ]; then
    echo "## Starting alignment with MAFFT..."
    set +u; conda activate /users/PAS0471/jelmer/.conda/envs/mafft-env; set -u
    
    mafft --reorder --auto --adjustdirection \
        --leavegappyregion $args_mafft "$fa_in" > "$fa_out"

elif [ "$align_program" = muscle ]; then
    echo "## Starting alignment with MUSCLE..."
    set +u; conda activate /users/PAS0471/jelmer/miniconda3/envs/muscle-env; set -u
    
    muscle -in "$fa_in" -out "$fa_out"

elif [ "$align_program" = tcoffee ]; then
    echo "## Starting alignment with t-coffee..."
    set +u; conda activate /fs/project/PAS0471/jelmer/conda/t_coffee-11.0.8; set -u

    fa_out_prefix=${fa_out%.*}
    
    t_coffee "$fa_in" -n_core="$n_cores" -output=fasta_aln -run_name="$fa_out_prefix"

    mv "$fa_out_prefix".fasta_aln "$fa_out"
else
    echo "## Not running alignment (alignment program variable set to: $align_program)"
fi

## Remove extra info after space from FASTA header lines
## ... and remove "_R_" prefixes for reverse-complemented seqs
sed -i -E -e 's/(^>[^ ]+) .*/\1/' -e 's/^>_R_/>/' "$fa_out"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n---------------------------"
echo -e "\n## Output files:"
ls -lh "$fa_out"
echo -e "\n## Done with script tree.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
