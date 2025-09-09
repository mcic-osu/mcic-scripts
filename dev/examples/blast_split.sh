#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --job-name=blast_split
#SBATCH --output=slurm-blast_split-%j.out

# Parse command-line args
fasta_glob=$1
blast_db=$2
outdir_base=$3
n_contigs_per_job=${4:-500}        # Optional: specify the number of contigs per job

# Load software
module load miniconda3
conda activate /fs/ess/PAS0471/jelmer/conda/seqkit
set -euo pipefail

# Report
echo "Starting script blast_split.sh"
date
echo "Input FASTA glob:         $fasta_glob"
echo "Number of input files:    $(ls $fasta_glob | wc -l)"
echo "BLAST db:                 $blast_db"
echo "Base output dir:          $outdir_base"
echo

# Create the output dir
mkdir -p "$outdir_base"/logs

# Split FASTA and run BLAST for each file
for input_fa in $fasta_glob; do
    nseqs=$(grep -c "^>" "$input_fa")
    echo -e "\n\n========Processing FASTA file $input_fa..."
    echo -e "Number of sequences: $nseqs"
    
    # Determine the output dirs
    file_id=$(basename "$input_fa" .fasta)
    outdir="$outdir_base"/"$file_id"
    fa_split_dir="$outdir"/split_fastas && mkdir -p "$fa_split_dir"
    
    # Split into files with N_CONTIGS_PER_JOB contigs each, for faster BLASTing
    echo -e "\nSplitting FASTA file..."
    seqkit split2 --force --by-size "$n_contigs_per_job" -O "$fa_split_dir" "$input_fa"

    # Loop over split FASTA files and run BLAST for each
    echo -e "\nSubmitting BLAST jobs..."
    for fa_split in "$fa_split_dir"/*part*.fasta; do
        file_id_split=$(basename "$fa_split" .fasta)
        outdir_split="$outdir"/"$file_id_split"
        
        # BLAST with evalue 1e-6 and only keeping the single best hit for each query (=input entry)
        sbatch -c8 mcic-scripts/blast/blast.sh \
            -i "$fa_split" \
            -o "$outdir_split" \
            --db "$blast_db" \
            --local \
            --top_n_query 1 \
            --evalue "1e-6"
    done
done

# Report
echo -e "\nDone with script blast_split.sh"
date
