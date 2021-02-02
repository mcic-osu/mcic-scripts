#!/bin/bash

# An example with primers 515F and 806R for 16S rRNA:

indir=data/fastq/raw
outdir=data/fastq/trimmed
primer_f=GAGTGYCAGCMGCCGCGGTAA
primer_r=ACGGACTACNVGGGTWTCTAAT

sbatch cutadapt.sh -i "$indir" -o "$outdir" -f "$primer_f" -r "$primer_r"
