#!/bin/bash
# Get reciprocal best BLAST hit between two sets of proteomes with DIAMOND

# Settings
g1_id=NCBI5
g2_id=Ensembl3
comp_id="$g1_id"_vs_"$g2_id"
PCT_ID=90
PCT_QCOV=80
PCT_SCOV=80

# Define input files
faa1_init=results/00_refdata/v5_NCBI/GCF_902167145.1.faa
gtf1=results/00_refdata/v5_NCBI/GCF_902167145.1.gtf
faa2_init=results/00_refdata/v3_ensembl/Zea_mays.AGPv3.22.pep.all.faa
gtf2=results/00_refdata/v3_ensembl/Zea_mays.AGPv3.22.gtf

# Define output files
outdir=results/genome_lookups/"$comp_id"
faa_dir="$outdir"/single-isoform-faa && mkdir -p "$faa_dir"
faa1="$faa_dir"/$(basename "$faa1_init")
faa2="$faa_dir"/$(basename "$faa2_init")
diamond_direction1="$outdir"/"$g1_id"_to_"$g2_id"/diamond_out.tsv
diamond_direction2="$outdir"/"$g2_id"_to_"$g1_id"/diamond_out.tsv
outfile="$outdir"/RBBH.tsv

# Create single-isoform protein files
sbatch mcic-scripts/annot/longest_tx.sh -i "$faa1_init" --gtf "$gtf1" -o "$faa_dir"
sbatch mcic-scripts/annot/longest_tx.sh -i "$faa2_init" --gtf "$gtf2" -o "$faa_dir" 

# Create DIAMOND databases
ls -lh "$faa1" "$faa2"
sbatch mcic-scripts/blast/diamond_db.sh -i "$faa1" -o "$outdir" --db_name "$g1_id".dmnd
sbatch mcic-scripts/blast/diamond_db.sh -i "$faa2" -o "$outdir" --db_name "$g2_id".dmnd

# Run DIAMOND
sbatch mcic-scripts/blast/diamond.sh \
    -i "$faa1" --db "$outdir"/"$g2_id".dmnd -o "$outdir"/"$g1_id"_to_"$g2_id" \
    --max_target_seqs 1 --pct_id "$PCT_ID" --pct_qcov "$PCT_QCOV" --pct_scov "$PCT_SCOV"
sbatch mcic-scripts/blast/diamond.sh \
    -i "$faa2" --db "$outdir"/"$g1_id".dmnd -o "$outdir"/"$g2_id"_to_"$g1_id" \
    --max_target_seqs 1 --pct_id "$PCT_ID" --pct_qcov "$PCT_QCOV" --pct_scov "$PCT_SCOV"

# Get reciprocal best blast hits
join -t $'\t' -1 1 -2 2 \
    <(tail -n+4 "$diamond_direction1" | cut -f 1,2 | sort -k1,1) \
    <(tail -n+4 "$diamond_direction2" | cut -f 1,2 | sort -k2,2) |
    awk '$2 == $3' |
    cut -f1,2  > "$outfile"
wc -l "$outfile"
