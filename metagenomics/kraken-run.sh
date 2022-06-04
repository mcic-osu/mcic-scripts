#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --mem=100G
#SBATCH --cpus-per-task=30
#SBATCH --output=slurm-kraken-run-%j.out

# HELP AND COMMAND-LINE OPTIONS ------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Run Kraken2 to assign taxonomy to sequences in a FASTA/FASTQ/pair of FASTQ file(s)"
    echo
    echo "## Syntax: $0 -i <input-file> -o <output-dir> -d <kraken-db-dir> ..."
    echo 
    echo "## Required options:"
    echo "  -i STRING           Input sequence file (FASTA, single-end FASTQ, or R1 from paired-end FASTQ)"
    echo "                          (If an R1 paired-end FASTQ file is provided, the name of the R2 file will be inferred.)"
    echo "  -o STRING           Output directory"
    echo "  -d STRING           Directory with an existing Kraken database"
    echo "                          (Use one of the scripts 'kraken-build-custom-db.sh' or 'kraken-build-std-db.sh' to create a Kraken database)"
    echo
    echo "## Other options:"
    echo "  -c NUM             Confidence required for assignment: number between 0 and 1          [default: 0.5]"
    echo "  -h                 Print this help message and exit"
    echo "  -m                 Don't load the full database into RAM memory                        [default: load into memory]"
    echo "                          (Can be useful for very large databases)"
    echo "  -n                 Add taxonomic names to the Kraken 'main' output file                [default: don't add]"
    echo "                          Note: this option is not compatible with Krona plotting)"
    echo "  -q INTEGER         Base quality Phred score required for use of a base in assignment   [default: 0]"
    echo "                          NOTE: If setting a score other than 0, any output sequence files (-w and -W options)"
    echo "                          will contain 'x's for masked bases."
    echo
    echo "  -w                 Write classified sequences to file                                  [default: don't write]"
    echo "  -W                 Write unclassified sequences to file                                [default: don't write]"
    echo
    echo "Example:          $0 -i data/A1_R1.fastq.gz -o results/kraken -d /fs/project/PAS0471/jelmer/refdata/kraken/PlusPFP"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
infile=""
outdir=""
krakendb_dir=""
add_names=false
min_conf=0.5
min_q=0
write_class=""
write_unclass=""
class_out_arg=""
unclass_out_arg=""
use_ram=true

## Get command-line options
while getopts 'i:o:d:c:q:mnwWh' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) krakendb_dir="$OPTARG" ;;
    c) min_conf="$OPTARG" ;;
    q) min_q="$OPTARG" ;;
    n) add_names=true ;;
    w) write_class=true ;;
    W) write_unclass=true ;;
    m) use_ram=false ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo "## Starting script kraken-run.sh..."
date
echo

## Process options
[[ "$infile" = "" ]] && echo "ERROR: must specify input file with -i" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "ERROR: input file $infile does note exist" >&2 && exit 1
[[ "$krakendb_dir" = "" ]] && echo "ERROR: must specify Kraken DB dir with -d" >&2 && exit 1
[[ ! -d "$krakendb_dir" ]] && echo "ERROR: Kraken DB dir $krakendb_dir does note exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "ERROR: must specify output dir with -o" >&2 && exit 1

[[ "$write_class" = true ]] && mkdir -p "$outdir"/classified
[[ "$write_unclass" = true ]] && mkdir -p "$outdir"/unclassified

# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env

## Bash strict settings
set -euo pipefail

## Report
echo "## Input file:                  $infile"
echo "## Output dir:                  $outdir"
echo "## Kraken db dir:               $krakendb_dir"
echo
echo "## Add tax. names:              $add_names"
echo "## Min. base qual:              $min_q"
echo "## Min. confidence:             $min_conf"
echo "## Use RAM memory to load database:    $use_ram"
echo

## Create output dir
mkdir -p "$outdir"

## Use RAM or not
if [[ "$use_ram" = false ]]; then
    mem_map_arg="--memory-mapping"
else
    mem_map_arg=""
fi

## Add tax. names or not -- when adding names, can't use the output for Krona plotting  
if [ "$add_names" = true ]; then
    names_arg="--use-names "
else
    names_arg=""
fi

## Make sure input file argument is correct based on file type 
if [[ "$infile" =~ \.fastq.gz$ ]]; then

    R1_in="$infile"
    R1_suffix=$(basename "$R1_in" | sed -E 's/.*(_R?[12]).*\.fa?s?t?q\.gz/\1/')
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    R1_basename=$(basename "$R1_in" .fastq.gz)
    sample_ID=${R1_basename/"$R1_suffix"/}

    if [[ -f $R2_in ]]; then
        echo "## Input is:                  paired FASTQ files"
        echo "## Input FASTQ file - R1:     $R1_in"
        echo "## Input FASTQ file - R2:     $R2_in"
        infile_arg="--gzip-compressed --paired $R1_in $R2_in"
        [[ ! -f "$R2_in" ]] && echo "## ERROR: R2 file $R2_in does not exist" >&2 && exit 1
        [[ "$R1_in" = "$R2_in" ]] && echo "## ERROR: R1 file $R1_in is the same as R2 file $R2_in" >&2 && exit 1

        if [[ "$write_class" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID#.fastq "
        fi
        if [[ "$write_unclass" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID#.fastq "
        fi

    else
        echo "## Input is:                  single-end FASTQ file"
        
        sample_ID=$(basename "$R1_in" .fastq.gz)
        
        infile_arg="--gzip-compressed $R1_in"

        if [[ "$write_class" = true ]]; then
            class_out_arg="--classified-out $outdir/classified/$sample_ID.fastq "
        fi
        if [[ "$write_unclass" = true ]]; then
            unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fastq "
        fi

    fi

else
    echo -e "## Input is:                   FASTA file"
    
    infile_basename=$(basename "$infile")
    sample_ID=${infile_basename%%.*}
    
    infile_arg="$infile"

    if [[ "$write_class" = true ]]; then
        class_out_arg="--classified-out $outdir/classified/$sample_ID.fa "
    fi
    if [[ "$write_unclass" = true ]]; then
        unclass_out_arg="--unclassified-out $outdir/unclassified/$sample_ID.fa "
    fi
fi

## Define output text files
outfile_main="$outdir"/"$sample_ID"_main.txt
outfile_report="$outdir"/"$sample_ID"_report.txt

## Report
echo "## Input file arg:            $infile_arg"
echo "## Sample ID:                 $sample_ID"
echo "## Output file - main:        $outfile_main"
echo "## Output file - report:      $outfile_report"
[[ "$write_class" = true ]] && echo "## Writing classified sequences:         $class_out_arg"
[[ "$write_unclass" = true ]] && echo "## Writing unclassified sequences:     $unclass_out_arg"
echo -e "------------------------\n"


# RUN KRAKEN -------------------------------------------------------------------
echo "## Starting Kraken2 run..."
kraken2 ${names_arg}--threads "$SLURM_CPUS_ON_NODE" \
    ${mem_map_arg}--minimum-base-quality "$min_q" \
    --confidence "$min_conf" \
    --report-minimizer-data \
    ${unclass_out_arg}--db "$krakendb_dir" \
    ${class_out_arg}--report "$outfile_report" \
    ${infile_arg}>"$outfile_main"

#? report-minimizer-data: see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information


# RENAME AND ZIP FASTQ FILES -------------------------
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

# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outfile_main" "$outfile_report"
[[ "$write_class" = true ]] && ls -lh "$outdir"/classified/"$sample_ID"*
[[ "$write_unclass" = true ]] && ls -lh "$outdir"/unclassified/"$sample_ID"*

echo -e "\n## Done with script kraken-run.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
