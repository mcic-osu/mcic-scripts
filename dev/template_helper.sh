#!/bin/bash

R1=
infile=

# FASTQ filename parsing
file_ext=$(basename "$R1" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2=${R1/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")

# Other input filename parsing
R1_filename="$(basename "$R1")"
file_ext=${R1_filename%.*}

# Make paths absolute
infile=$(realpath "$infile")
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir" 

# Read a fofn
[[ -n "$fofn" ]] && mapfile -t infiles <"$fofn"
[[ -n "$indir" ]] && mapfile infiles < <(find "$indir" -type f)

# Check derived inputs
[[ ! -f $R2 ]] && die "Input file R2 ($R2) does not exist"
[[ "$R1" == "$R2" ]] && die "Input R1 and R2 FASTQ files are the same file: $R1"
