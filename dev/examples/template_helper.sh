#!/bin/bash

R1=
infile=

# FASTQ filename parsing
file_ext=$(basename "$R1" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E 's/.*(_R?[12]).*/\1/')
R2_suffix=${R1_suffix/1/2}
R2=$(echo "$R1" | sed -E "s/${R1_suffix}([._])/${R2_suffix}\1/")
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")

# Other input filename parsing
R1_filename="$(basename "$R1")"
file_ext=${R1_filename##*.}

# Make paths absolute
infile=$(realpath "$infile")
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir" 

# Read a fofn
[[ -n "$fofn" ]] && mapfile -t infiles <"$fofn"

# Create an array of input files
[[ -n "$indir" ]] && mapfile infiles < <(find "$indir" -type f)
mapfile -t fastas < <(find "$indir" -iname '*.fasta' -or -iname '*.fa' -or -iname '*.fna' -or -iname '*.fna.gz')
[[ ${#fastas[@]} -eq 0 ]] && die "No FASTA files found..."

# Check derived inputs
[[ "$R1" == "$R2" ]] && die "Input R1 and R2 FASTQ files are the same file: $R1"

# Process positional args as input files (for an example, see misc/fastqc.sh)
declare -a infiles
count=0
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && readonly infiles[0]=$1 ;;
        * )                 infiles[$count]=$1 && count=$(( count + 1 )) ;;
    esac
    shift
done
[[ "${#infiles[@]}" -eq 0 ]] && die "Please specify input file(s) with -i/--infile or as positional arguments" "$all_opts"
for infile in "${infiles[@]}"; do
    [[ ! -f "$infile" ]] && die "Input file $infile does not exist"
done
echo "Number of input files:            ${#infiles[@]}"
ls -lh "${infiles[@]}"
