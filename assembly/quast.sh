#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --job-name=quast
#SBATCH --output=slurm-quast-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "         Run QUAST to check the quality of a genome assembly"
    echo "======================================================================"
    echo "USAGE:"
    echo "sbatch $0 [ -i <input-file> / -d <input-dir> ] -o <output-dir> ..."
    echo
    echo "REQUIRED OPTIONS:"
    echo "  Specify the input with EITHER '-i' (single file) or '-d' (all files in a dir)"
    echo "  -i <file>         Input file (genome assembly nucleotide FASTA)"
    echo "  -d <dir>          Input dir (Only files with the extension '.fasta' will be included)"
    echo "  -o <dir>          Output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -r <file>         Reference genome nucleotide FASTA file    [default: no reference genome]"
    echo "  -g <file>         Reference genome GFF file                 [default: no reference genome]"
    echo "  -1 <file>         FASTQ file with forward reads             [default: no reads]"
    echo "  -2 <file>         FASTQ file with reverse reads             [default: no reads]"
    echo "  -a <string>       Quoted string with additional argument(s) to pass to QUAST"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h              Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly/my.fasta -o results/quast"
    echo "  sbatch $0 -d results/assemblies -o results/quast"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - GitHub repo: https://github.com/ablab/quast"
    echo "  - Manual: https://quast.sourceforge.net/docs/manual.html"
    echo "  - Paper: https://academic.oup.com/bioinformatics/article/29/8/1072/228832"
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    source activate /fs/ess/PAS0471/jelmer/conda/quast-5.0.2
}

## Print version
Print_version() {
    Load_software
    quast.py --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
dryrun=false
debug=false

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
indir=""
outdir=""
ref_fna="" && ref_fna_arg=""
ref_gff="" && ref_gff_arg=""
R1="" && R1_arg=""
R2="" && R2_arg=""
more_args=""

## Parse command-line options
while getopts '1:2:i:d:o:r:g:a:h' flag; do
    case "${flag}" in
        i) infile="$OPTARG" ;;
        d) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        1) R1="$OPTARG" ;;
        2) R2="$OPTARG" ;;
        r) ref_fna="$OPTARG" ;;
        g) ref_gff="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Print_help; exit 0 ;;
        \?) Die "Invalid option -$OPTARG" ;;
        :) Die "Option -$OPTARG requires an argument" ;;
    esac
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## Load the software
Load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ "$indir" = "" && "$infile" = "" ]] && Die "Please specify either an input file with -i, or an input dir with -d"
[[ "$indir" != "" && "$infile" != "" ]] && Die "Please specify either an input file with -i, or an input dir with -d, and not both!"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o"
[[ "$R2" != "" && "$R1" = "" ]] && Die "When specifying a FASTQ file with reverse reads ('-2'), please also specify a FASTQ file with forward reads ('-1')"
[[ "$infile" != "" && ! -f "$infile" ]] && Die "Input file (-i) $infile does not exist or is not a regular file"
[[ "$indir" != "" && ! -d "$indir" ]] && Die "Input dir (-d) $indir does not exist or is not a directory"

## Build argument for assembly input
if [[ $indir != "" ]]; then
    infile_arg="$indir/*.fasta"
else
    infile_arg="$infile"
fi

## If the input is a single file, make a separate output dir
if [[ $infile != "" ]]; then
    sampleID=$(basename "$infile" | sed -E 's/.fn?as?t?a?//')
    outdir=$outdir/"$sampleID"
fi

## Build arguments
[[ $R1 != "" ]] && R1_arg="-1 $R1"
[[ $R2 != "" ]] && R2_arg="-2 $R2"
[[ $ref_fna != "" ]] && ref_fna_arg="-r $ref_fna"
[[ $ref_gff != "" ]] && ref_gff_arg="-g $ref_gff"

## Get number of threads
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT QUAST.SH"
date
echo "=========================================================================="
[[ $infile != "" ]] && echo "Input file:                           $infile"
[[ $indir != "" ]] && echo "Input dir:                            $indir"
echo "Output dir:                           $outdir"
[[ $ref_fna != "" ]] && echo "Reference FASTA file:                 $ref_fna"
[[ $ref_gff != "" ]] && echo "Reference GFF file:                   $ref_gff"
[[ $R1 != "" ]] && echo "R1 FASTQ file:                        $R1"
[[ $R2 != "" ]] && echo "R2 FASTQ file:                        $R2"
[[ $more_args != "" ]] && echo "Other arguments to pass to QUAST:     $more_args"
if [[ $indir != "" ]]; then
    echo -e "\nListing input FASTA files:"
    ls -lh "$indir"/*.fasta
fi
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
## Make output dir
mkdir -p "$outdir"

[[ "$debug" = false ]] && set -o xtrace

quast.py \
    --conserved-genes-finding \
    -o "$outdir" \
    -t "$n_threads" \
    $ref_fna_arg \
    $ref_gff_arg \
    $R1_arg \
    $R2_arg \
    $more_args \
    $infile_arg

[[ "$debug" = false ]] && set +o xtrace


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing files in the output dir:"
    ls -lh "$outdir"
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
echo "## Done with script"
date
