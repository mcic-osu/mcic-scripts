#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --job-name=fmlrc2
#SBATCH --output=slurm-fmlrc2-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                  $0"
    echo " RUN FMLRC2 TO CORRECT LONG READS OR AN ASSEMBLY WITH ILLUMINA READS"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <assembly/long reads> -I <short-read FASTQ> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input assembly (FASTA) or long-read gzipped FASTQ file"
    echo "  -I/--fq_short   <str>   Input short-read FASTQ file(s)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --eukaryote             Turn on 'eukaryote mode'                    [default: off]"
    echo "                          See https://www.biorxiv.org/content/10.1101/2022.07.22.501182v1 for details"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to FMLRC2"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v/--version            Print the version of FMLRC2 and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly.fasta -I data/my.fastq -o results/fmlrc2 "
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - ..."
    echo
    echo "OUTPUT:"
    echo "  - ..."
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - Repo/docs: https://github.com/HudsonAlpha/fmlrc2"
    echo "  - Paper: https://www.biorxiv.org/content/10.1101/2022.07.22.501182v1.abstract"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    source activate /fs/ess/PAS0471/jelmer/conda/fmlrc2-0.1.7 # Also has seqtk installed
}

## Print version
Print_version() {
    Load_software
    fmlrc2 --version
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
## Option defaults
eukaryote=false && eukar_arg=""

debug=false
dryrun=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
infile=""
fq_short=""
outdir=""
more_args=""

## Parse command-line args
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )         shift && infile=$1 ;;
        -I | --fq_short )       shift && fq_short=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --eukaryote )           eukaryote=true ;;
        --more_args )           shift && more_args=$1 ;;
        --debug )               debug=true ;;
        --dryrun )              dryrun=true ;;
        -v | -v | --version )        Print_version; exit ;;
        -h | --help )           Print_help; exit ;;
        * )                     Print_help; Die "Invalid option $1" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## Get number of threads
if [[ "$dryrun" = false ]]; then
    if [[ -z "$SLURM_CPUS_PER_TASK" ]]; then
        n_threads="$SLURM_NTASKS"
    else
        n_threads="$SLURM_CPUS_PER_TASK"
    fi
fi

## Bash script settings
set -euo pipefail

## Check input
[[ $infile = "" ]] && Die "Please specify an assembly or long-read FASTQ file with -i/--infile"
[[ $fq_short = "" ]] && Die "Please specify one more short-read FASTQ files with -I/--fq_short"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ ! -f $infile ]] && Die "Input file $infile does not exist"

## Define output files
bwt_out="$outdir"/comp_msbwt.npy

file_ext=$(echo "$infile" | sed -E 's/.*(\.fasta|\.fastq.gz|\.fq\.gz)/\1/')
file_id=$(basename "$infile" "$file_ext")
fasta_out="$outdir"/"$file_id".fasta

if [[ "$infile" =~ q.gz$ ]]; then
    ## FASTQ
    init_command="gunzip -c $infile | awk 'NR % 4 == 2'"
else
    ## FASTA (Linearize with 'seqtk seq')
    init_command="cat $infile | seqtk seq - | grep -v '^>'"
fi

[[ "$eukaryote" = true ]] && eukar_arg="-k 21 -k 59 -k 80 --min_dynamic_count 0"

## Report
echo
echo "=========================================================================="
echo "               STARTING SCRIPT FMLRC2.SH"
date
echo "=========================================================================="
echo "Input assembly or long-read FASTQ file:   $infile"
echo "Input short-read FASTQ file(s):           $fq_short"
echo "Output dir:                               $outdir"
echo "Eukaryote mode:                           $eukaryote"
echo
echo "Output BWT file:                          $bwt_out"
echo "Output FASTA file:                        $fasta_out"
[[ $more_args != "" ]] && echo "Other arguments for FMLRC2:    $more_args"
[[ $dryrun = true ]] && echo "THIS IS A DRY-RUN"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    ## Create the output directory
    mkdir -p "$outdir"/logs

    [[ "$debug" = false ]] && set -o xtrace

    echo -e "\n## Now building the BWT with ropebwt2..."
    if [[ ! -f "$bwt_out" ]]; then
        eval $init_command |
            tr NT TN |
            ropebwt2 -LR |
            tr NT TN |
            fmlrc2-convert "$bwt_out"
    else
        echo -e "\nFile $bwt_out detected, not rerunning fmlrc2-convert!"
        ls -lh "$bwt_out"
    fi

    echo -e "\n## Now running fmlrc2..."
    date
    fmlrc2 \
        --threads "$n_threads" \
        --cache_size 11 \
        $eukar_arg \
        $more_args \
        "$bwt_out" \
        "$infile" \
        "$fasta_out"

    [[ "$debug" = false ]] && set +o xtrace
fi

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
