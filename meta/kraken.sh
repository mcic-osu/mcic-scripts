#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --mem=100G
#SBATCH --cpus-per-task=30
#SBATCH --job-name=kraken
#SBATCH --output=slurm-kraken-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "===================================================================================="
    echo "                  $0"
    echo "Run Kraken2 to assign taxonomy to sequences in a FASTA/FASTQ/pair of FASTQ file(s)"
    echo "===================================================================================="
    echo
    echo "USAGE:"
    echo "sbatch $0 -i <input-file> -o <output-dir> -d <kraken-db-dir> ..."
    echo "bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i FILE            Input sequence file (FASTA, single-end FASTQ, or R1 from paired-end FASTQ)"
    echo "                       - If an R1 paired-end FASTQ file is provided, the name of the R2 file will be inferred"
    echo "                       - FASTA files should be unzipped; FASTQ files should be gzipped"
    echo "  -o DIR             Output directory"
    echo "  -d DIR             Directory with an existing Kraken database"
    echo "                          (Use one of the scripts 'kraken-build-custom-db.sh' or 'kraken-build-std-db.sh' to create a Kraken database)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -c NUM             Confidence required for assignment: number between 0 and 1            [default: 0.5]"
    echo "  -q INTEGER         Base quality Phred score required for use of a base in assignment     [default: 0]"
    echo "                          NOTE: If setting a score other than 0, any output sequence files (-w and -W options)"
    echo "                                will contain 'x's for masked bases."
    echo "  -m                 Don't load the full database into RAM memory                          [default: load into memory]"
    echo "                          (Considerably lower, but can be useful/needed with very large databases)"
    echo "  -n                 Add taxonomic names to the Kraken 'main' output file                  [default: don't add]"
    echo "                          NOTE: This option is not compatible with Krona plotting"
    echo "  -s                 FASTQ files are single-end                                            [default: paired-end]"
    echo "  -w                 Write 'classified' reads/sequences to file (in '<outdir>/classified' dir)     [default: don't write]"
    echo "  -W                 Write 'unclassified' reads/sequences to file (in '<outdir>/unclassified' dir) [default: don't write]"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -x               Run the script in debug mode"
    echo "    -h               Print this help message and exit"
    echo "    -v               Print the version of Kraken and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "sbatch $0 -i data/A1_R1.fastq.gz -o results/kraken -d /fs/project/PAS0471/jelmer/refdata/kraken/std"
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env
}

## Print version
Print_version() {
    Load_software
    kraken2 --version
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
min_conf=0.5
min_q=0
add_names=false && names_arg=""
use_ram=true && mem_map_arg=""
single_end=false

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder variables
infile=""
outdir=""
krakendb_dir=""
write_class="" && class_out_arg=""
write_unclass="" && unclass_out_arg=""
dryrun=false
debug=false

## Get command-line options
while getopts 'i:o:d:c:q:sNXmnwWh' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) krakendb_dir="$OPTARG" ;;
    c) min_conf="$OPTARG" ;;
    q) min_q="$OPTARG" ;;
    n) add_names=true ;;
    s) single_end=true ;;
    w) write_class=true ;;
    W) write_unclass=true ;;
    m) use_ram=false ;;
    N) dryrun=true ;;
    X) debug=true ;;
    h) Print_help && exit 0 ;;
    \?) Die "Invalid option -$OPTARG" ;;
    :) Die "Option -$OPTARG requires an argument" ;;
    esac
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
[[ "$debug" = true ]] && set -o xtrace

## Load software
[[ "$dryrun" = false ]] && Load_software

## Bash strict settings
set -euo pipefail

## Report
echo
echo "=========================================================================="
echo "                     STARTING SCRIPT KRAKEN.SH"
date
echo "=========================================================================="
echo "Input file:                     $infile"
echo "Output dir:                     $outdir"
echo "Kraken db dir:                  $krakendb_dir"
echo
echo "Add tax. names:                 $add_names"
echo "Min. base qual:                 $min_q"
echo "Min. confidence:                $min_conf"
echo "Use RAM to load database:       $use_ram"
echo

## Process options
[[ "$infile" = "" ]] && Die "Must specify input file with -i"
[[ ! -f "$infile" ]] && Die "Input file $infile does note exist"
[[ "$krakendb_dir" = "" ]] && Die "Must specify a Kraken DB dir with -d"
[[ ! -d "$krakendb_dir" ]] && Die "Kraken DB dir $krakendb_dir does note exist"
[[ "$outdir" = "" ]] && Die "Must specify an output dir with -o"
[[ "$min_conf" = "" ]] && Die "Min confidence is not set"
[[ "$min_q" = "" ]] && Die "Min quality is not set"

[[ "$write_class" = true ]] && mkdir -p "$outdir"/classified
[[ "$write_unclass" = true ]] && mkdir -p "$outdir"/unclassified

## Build Kraken args (leave space after!)
### RAM
[[ "$use_ram" = false ]] && mem_map_arg="--memory-mapping "
### Add tax. names or not -- when adding names, can't use the output for Krona plotting  
[[ "$add_names" = true ]] && names_arg="--use-names "

## Make sure input file argument is correct based on file type 
if [[ "$infile" =~ \.fa?s?t?q.gz$ ]]; then

    R1_in="$infile"
    R1_basename=$(basename "$R1_in" | sed -E 's/\.fa?s?t?q\.gz//')
    R1_suffix=$(echo "$R1_basename" | sed -E 's/.*(_R?[12]).*/\1/')
    
    if [[ "$single_end" = false ]]; then

        echo "Input type is:                  paired-end FASTQ files"
        
        if [[ "$R1_suffix" != "$R1_basename" && "$R1_suffix" != "" ]]; then
            die "Can't figure out R2 filename"
        fi

        R2_suffix=${R1_suffix/1/2}
        R2_in=${R1_in/$R1_suffix/$R2_suffix}
        sample_ID=${R1_basename/"$R1_suffix"/}
        infile_arg="--gzip-compressed --paired $R1_in $R2_in"

        echo "Input FASTQ file - R1:          $R1_in"
        echo "Input FASTQ file - R2:          $R2_in"

        [[ ! -f "$R2_in" ]] && Die "R2 file $R2_in does not exist"
        [[ "$R1_in" = "$R2_in" ]] && Die "R1 file $R1_in is the same as R2 file $R2_in"

        if [[ "$write_class" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID#.fastq "
        fi
        if [[ "$write_unclass" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID#.fastq "
        fi

    else

        echo "Input type is:                  single-end FASTQ file"
        sample_ID=$(basename "$R1_in" .fastq.gz)
        infile_arg="--gzip-compressed $R1_in"

        if [[ "$write_class" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID.fastq "
        fi
        if [[ "$write_unclass" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fastq "
        fi
    
    fi

elif [[ "$infile" =~ \.fn?a?s?t?a$ ]]; then
    echo -e "Input type is:                   FASTA file"
    infile_basename=$(basename "$infile")
    sample_ID=${infile_basename%%.*}
    infile_arg="$infile"

    if [[ "$write_class" = true ]]; then
        class_out_arg="--classified-out $outdir/classified/$sample_ID.fa "
    fi
    if [[ "$write_unclass" = true ]]; then
        unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fa "
    fi

else
    Die "Unknown input file type"
fi

## Define output text files
outfile_main="$outdir"/"$sample_ID"_main.txt
outfile_report="$outdir"/"$sample_ID"_report.txt

## Report
echo "Sample ID (inferred):           $sample_ID"
echo "Output file - main:             $outfile_main"
echo "Output file - report:           $outfile_report"
[[ "$write_class" = true ]] && echo "Writing classified sequences:   $class_out_arg"
[[ "$write_unclass" = true ]] && echo "Writing unclassified sequences: $unclass_out_arg"
echo "=========================================================================="
echo


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    ## Create output dir
    mkdir -p "$outdir"/logs

    ## Run Kraken
    echo "## Starting Kraken2 run..."
    [[ "$debug" = false ]] && set -o xtrace

    kraken2 ${names_arg}--threads "$SLURM_CPUS_ON_NODE" \
        ${mem_map_arg}--minimum-base-quality "$min_q" \
        --confidence "$min_conf" \
        --report-minimizer-data \
        ${unclass_out_arg}--db "$krakendb_dir" \
        ${class_out_arg}--report "$outfile_report" \
        ${infile_arg}> "$outfile_main"

    [[ "$debug" = false ]] && set +o xtrace

    #? report-minimizer-data: see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information

    ## Rename and zip FASTQ files
    if [[ "$write_class" = true ]]; then
        mv "$outdir"/classified/"$sample_ID"_1.fastq "$outdir"/classified/"$sample_ID"_R1.fastq
        gzip -f "$outdir"/classified/"$sample_ID"_R1.fastq

        mv "$outdir"/classified/"$sample_ID"_2.fastq "$outdir"/classified/"$sample_ID"_R2.fastq
        gzip -f "$outdir"/classified/"$sample_ID"_R2.fastq
    fi

    if [[ "$write_unclass" = true ]]; then
        mv "$outdir"/unclassified/"$sample_ID"_1.fastq "$outdir"/unclassified/"$sample_ID"_R1.fastq
        gzip -f "$outdir"/unclassified/"$sample_ID"_R1.fastq

        mv "$outdir"/unclassified/"$sample_ID"_2.fastq "$outdir"/unclassified/"$sample_ID"_R2.fastq
        gzip -f "$outdir"/unclassified/"$sample_ID"_R2.fastq
    fi
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "## Done with script"
date
if [[ "$dryrun" = false ]]; then
    echo -e "\n## Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n## Listing output files:"
    ls -lh "$outfile_main" "$outfile_report"
    [[ "$write_class" = true ]] && ls -lh "$outdir"/classified/"$sample_ID"*
    [[ "$write_unclass" = true ]] && ls -lh "$outdir"/unclassified/"$sample_ID"*
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
fi
echo
