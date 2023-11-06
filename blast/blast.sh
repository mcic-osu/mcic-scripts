#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --job-name=blast
#SBATCH --output=slurm-blast-%j.out

#TODO - Add taxonomic info with taxonkit

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# BLAST docs: https://www.ncbi.nlm.nih.gov/blast/BLAST_guide.pdf

# Constants - generic
DESCRIPTION="Run NCBI BLAST on an input (query) FASTA file, and optionally download
aligned sequences and/or genomes. The input (query) FASTA file can contain multiple
or even many sequences, though it will be quicker to split a multiFASTA file,
and submit a separate job for each single-sequence FASTA file.
Additionally, downloaded sequences (e.g. with --download_genomes) are currently not
output separately for each query.

OUTPUT:
  All output will be placed inside the specified output dir.
  The output will include TSV files with raw ('blast_out_raw.tsv') and filtered
  ('blast_out_filtered.tsv'), and if requested, downloaded sequences in separate
  subdirectories."
SCRIPT_VERSION="2023-10-29"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
VERSION_COMMAND="blastn -version; datasets --version"
export TOOL_NAME="NCBI BLAST+ and datasets" 
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
export LC_ALL=C                     # Locale for sorting

# Constants - settings
# - With genome download, get separate metadata file from the NCBI datasets tool with the following fields:
META_FIELDS="accession,assminfo-name,organism-name,assminfo-refseq-category,assminfo-level,assmstats-number-of-contigs,assmstats-contig-n50"
# - Am not able to get scientific name of subject seq to be included ('ssciname' / 'sscinames')
# - In addition to the 'qcovhsp' included above, there is also 'qcovs', which will contain the total coverage across all HSPs
DEFAULT_LOCAL_DB_NT=/fs/project/pub_data/blast-database/2023-06/nt              # Default local DB for nucleotide BLAST
DEFAULT_LOCAL_DB_AA=/fs/project/pub_data/blast-database/2023-06/nr              # Default local DB for protein BLAST

# Defaults - generic
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/blast
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                  # When true, just print tool & script version info and exit

# Defaults - settings
blast_format="6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp slen stitle staxids"
local=false && remote_opt=" -remote" # Run BLAST locally (own db) or remotely (NCBI's db over the internet)
remote_db_nt=nt                     # Default remote db for BLAST-to-nucleotide (blastn, tblastx, tblastn)
remote_db_aa=nr                     # Default remote db for BLAST-to-protein (blastp and blastx)
blast_type=blastn                   # BLAST type
db_type=nuc                         # 'prot' (proteins/amino acids) or 'nuc' (nucleotides)  (automatically determined)
dl_db=nuccore                       # 'nuccore' for nucleotide db, 'protein' for protein db (automatically determined)
blast_task=                         # 'task' within BLAST type, e.g. 'megablast' for blastn
top_n_query=100                     # Keep the top-N hits only for each query (empty => keep all)
top_n_subject=                      # Keep the top-N hits only for each subject, per query (empty => keep all)
evalue="1e-6"                       # E-value threshold
pct_id=                             # % identity threshold (empty => no threshold)
pct_cov=                            # Threshold for % of query covered by the alignment length (empty => no threshold)
force=false                         # Don't rerun BLAST if the output file already exists
to_find_genomes=false               # Find genomes of subjects and create lookup tables with subject acessions?
to_dl_genomes=false                 # Download full genomes of subjects?
to_dl_subjects=false                # Download full subjects?
to_dl_aligned=false                 # Download aligned sequences?
add_header=true                     # Add column header to final BLAST output file

# ==============================================================================
#                           GENERIC FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage - will run BLAST remotely with the nt database & no output filtering or sequence downloading:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast"
    echo
    echo "  - Run a local BLAST using OSCs June 2023 BLAST db (used by default with --local):"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --local"
    echo
    echo "  - Limit online BLAST database to specific taxa (using NCBI taxon IDs):"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --tax_ids '343,56448'"
    echo
    echo "  - Download aligned parts of sequences, full subjects, and full genomes:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --dl_aligned --dl_subjects --dl_genomes"
    echo
    echo "  - Use % identity and query coverage thresholds:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --pct_id 90 --pct_cov 90"
    echo
    echo "  - Keep only the best 10 hits per query:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --top_n 10"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input FASTA file (can contain one or more sequences)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "GENERAL OPTIONS (OPTIONAL):"
    echo "  --no_header                 Don't add column headers to final BLAST output TSV file [default: add]"
    echo "                                The header won't be added to the raw output file, which can be used for filtering"
    echo "  --force                     Run BLAST even if the output file already exists        [default: $force]"
    echo "                                When false, BLAST won't be rerun but downstream operations will be"
    echo
    echo "GENERAL BLAST OPTIONS (OPTIONAL):"
    echo "  --local                     Run BLAST with a local (on-disk) database               [default: $local]"
    echo "  --db                <str>   - If running remotely: NCBI database name like 'nt'/'nr'[default: 'nt' for nucleotide, 'nr' for protein]"
    echo "                              - If running locally: default is 'nt' or 'nr' from '/fs/project/pub_data/blast-database/2023-06'"
    echo "                              - To run with a custom local db, use '--local' AND specify the prefix (dir + db name, no file extensions) of a local BLAST db"
    echo "  --blast_type        <str>   BLAST type: 'blastn', 'blastp', 'blastx', 'tblastx', or 'tblastn' [default: $blast_type]"
    echo "  --blast_task        <str>   'Task' for blastn or blastp, e.g. 'megablast' for blastn  [default: BLAST program default]"
    echo "                              For blastn, the default is 'megablast', other options are: 'blastn', 'blastn-short', 'dc-megablast', 'rmblastn'"
    echo "                              For blastp, the default is 'blastp', other options are: 'blastp-fast', 'blastp-short'"
    echo "                                See https://www.ncbi.nlm.nih.gov/books/NBK569839/#usrman_BLAST_feat.Tasks"
    echo "  --out_format        <str>   BLAST output format string. NOTE: changing this may mess up output filtering steps, which rely on the default format"
    echo "                                [default: $blast_format]"
    echo
    echo "BLAST THRESHOLD AND FILTERING OPTIONS (OPTIONAL):"
    echo "  --tax_ids           <str>   Comma-separated list of NCBI taxon IDs (just the numbers, no 'txid' prefix)"
    echo "                                The BLAST search will be limited to these taxa        [default: use full database]"
    echo "                                NOTE: This option only works for remote, nucleotide-based searches!"
    echo "  --max_target_seqs   <int>   Max. nr of target sequences to keep                     [default: BLAST default (=500)]"
    echo "                                This option will be applied *during* the BLAST run"
    echo "                                This number should be increased from the default when ~hundreds of hits are expected"
    echo "  --evalue            <num>   E-value threshold in scientific notation                [default: $evalue]"
    echo "                                This option will be applied *during* the BLAST run"
    echo "  --pct_id            <num>   Percentage identity threshold                           [default: none]"
    echo "                                This threshold will be applied *after* running BLAST"
    echo "  --pct_cov           <num>   Threshold for % of query covered by the alignment       [default: none]"
    echo "                                This threshold will be applied *after* running BLAST"
    echo "  --top_n_query       <int>   Only keep the top N hits for each query                 [default: $top_n_query]"
    echo "                                This threshold will be applied *after* running BLAST"
    echo "                                A threshold of 0 means no filtering"
    echo "  --top_n_subject     <int>   Only keep the top N hits for each subject, per query    [default: keep all]"
    echo "                                This threshold will be applied *after* running BLAST"
    echo "                                A threshold of 0 means no filtering"
    echo
    echo "SEQUENCE LOOKUP AND DOWNLOAD OPTIONS (OPTIONAL):"
    echo "  --dl_aligned                Download aligned parts of subject (db) sequences        [default: $to_dl_aligned]"
    echo "  --dl_subjects               Download full subject (db) sequences that were aligned  [default: $to_dl_subjects]"
    echo "                                For protein BLAST, this will also download nucleotide sequences for non 'WP_' (multi-species) accessions"
    echo "  --find_genomes              Find genome accession numbers of aligned sequences      [default: $to_find_genomes]"
    echo "  --dl_genomes                Download full genomes of aligned sequences              [default: $to_dl_genomes]"
    echo
    echo "UTILITY OPTIONS (OPTIONAL):"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you must provide one to run the script with a container.)"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --no_strict                 Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of BLAST and the NCBI datasets tool and exit"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    # Determine the location of this script, and based on that, the function script
    if [[ "$IS_SLURM" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/"$(basename "$FUNCTION_SCRIPT_URL")")
    # Download the function script if needed, then source it
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        function_script=$(basename "$FUNCTION_SCRIPT_URL")
        wget "$FUNCTION_SCRIPT_URL" -O "$function_script"
    fi
    source "$function_script"
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                           ANALYSIS FUNCTIONS
# ==============================================================================
run_blast() {
    log_time "Now running BLAST..."
    #(Options need to be awkwardly collapsed like this or BLAST will choke on the empty spaces)
    runstats "$blast_type" \
        -db "$db" \
        -query "$infile" \
        -out "$blast_out_raw" \
        -outfmt "$blast_format" \
        -evalue "$evalue" \
        ${maxtarget_opt}${task_opt}${remote_opt}${thread_opt}${tax_opt}${spacer}"${tax_optarg}"
}

process_blast() {
    log_time "Nr of hits in the raw BLAST output file: $(wc -l < "$blast_out_raw")"

    # Sort by: query name (1), then e-value (11), then bitscore (12), then % identical (3)
    log_time "Sorting BLAST output by goodness of the match"
    sort -k1,1 -k11,11g -k12,12gr -k3,3gr "$blast_out_raw" > "$blast_out_sorted"

    # 1. Filter by percent identity
    if [[ -n "$pct_id" ]]; then
        log_time "Filtering output using a percent identity threshold of $pct_id"
        blast_out_id="$outdir"/blast_out_pctid.tsv
        
        awk -F"\t" -v OFS="\t" -v pct_id="$pct_id" \
            '$3 >= pct_id' "$blast_out_sorted" > "$blast_out_id"
        log_time "Retained $(wc -l < "$blast_out_id") of $(wc -l < "$blast_out_sorted") hits"
    else
        blast_out_id="$blast_out_sorted"
    fi

    # 2. Filter by percent coverage
    if [[ -n "$pct_cov" ]]; then
        log_time "Filtering output using a percent coverage threshold of $pct_cov"
        blast_out_cov="$outdir"/blast_out_cov.tsv

        awk -F"\t" -v OFS="\t" -v pct_cov="$pct_cov" \
            '$14 >= pct_cov' "$blast_out_id" > "$blast_out_cov"
        log_time "Retained $(wc -l < "$blast_out_cov") of $(wc -l < "$blast_out_id") hits"
    else
        blast_out_cov="$blast_out_id"
    fi

    # 3. Only retain top-N matches per query
    if [[ -n "$top_n_query" && "$top_n_query" != 0 ]]; then
        log_time "Getting the top $top_n_query hits for each query"
        blast_out_topq="$outdir"/blast_out_topq.tsv

        while read -r query; do
            grep -w -m "$top_n_query" "$query" "$blast_out_cov"
        done < <(cut -f1 "$blast_out_cov" | sort -u) |
            sort -k1,1 -k11,11g -k12,12gr -k3,3gr > "$blast_out_topq"
        
        log_time "Retained $(wc -l < "$blast_out_topq") of $(wc -l < "$blast_out_cov") hits"
    else
        blast_out_topq="$blast_out_cov"
    fi

    # 4. Only retain top-N matches per subject
    if [[ -n "$top_n_subject" && "$top_n_subject" != 0 ]]; then
        log_time "Getting the top $top_n_subject hits for each subject"
        
        while read -r subject; do
            grep -w -m "$top_n_subject" "$subject" "$blast_out_topq"
        done < <(cut -f2 "$blast_out_topq" | sort -u) |
        sort -k1,1 -k11,11g -k12,12gr -k3,3gr > "$blast_out_final"
        
        log_time "Retained $(wc -l < "$blast_out_final") of $(wc -l < "$blast_out_topq") hits"
    else
        cp "$blast_out_topq" "$blast_out_final"
    fi

    # Clean & report
    [[ -f "$blast_out_cov" ]] && rm "$blast_out_cov"
    [[ -f "$blast_out_id" ]] && rm "$blast_out_id"
    [[ -f "$blast_out_sorted" ]] && rm "$blast_out_sorted"
    [[ -f "$blast_out_topq" ]] && rm "$blast_out_topq"

    log_time "Listing the final BLAST output file:"
    ls -lh "$blast_out_final"

    log_time "Showing the first few lines of the final BLAST output file (without header):"
    head -n 5 "$blast_out_final"

    n_subjects=$(cut -f 2 "$blast_out_final" | sort -u | wc -l)
    n_queries=$(cut -f 1 "$blast_out_final" | sort -u | wc -l)
    log_time "Number of distinct subjects in the final BLAST output file: $n_subjects"
    log_time "Number of distinct queries in the final BLAST output file: $n_queries"
}

find_genomes() {
    echo -e "\n================================================================"
    log_time "Now finding full genomes of matched sequences..."
    mkdir -p "$outdir"/genomes

    # Move into outdir & define output files
    assembly_list="$outdir"/genomes/assemblies.txt      # For 'datasets' download command
    assembly_lookup="$outdir"/genomes/assemblies.tsv
    meta_file="$outdir"/genomes/genome_metadata.tsv

    # Get the RefSeq (GCF_) assembly ID
    log_time "Looking up the genome accession IDs..."
    > "$assembly_list"
    > "$assembly_lookup"
    while read -r -u 9 accession; do
        # Temporary accession list file
        assembly_list_acc="$outdir"/genomes/tmp_"$accession".txt

        if [[ "$dl_db" == "nuccore" ]]; then
            # For nucleotide searches

            # First attempt to get assembly ID
            #(Note: 'head' at end because in some cases, multiple assembly versions are associated with an accession)
            mapfile -t assemblies < <(esearch -db "$dl_db" -query "$accession" | elink -target assembly |
                esummary | xtract -pattern DocumentSummary -element RefSeq | head -n 1)
            
            # If needed, second attempt to get assembly ID
            if [[ ${#assemblies[@]} -eq 0 ]]; then
                mapfile -t assemblies < <(esearch -db "$dl_db" -query "$accession" | elink -target assembly |
                    esummary | xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)
            fi
        else
            # For protein searches
            if [[ ! "$accession" =~ ^WP_ ]]; then
                # Non-'WP_' (multispecies) entries
                mapfile -t assemblies < <(esearch -db "$dl_db" -query "$accession" |
                    elink -target nuccore | elink -target assembly | esummary |
                    xtract -pattern DocumentSummary -element RefSeq | head -n 1)
                
                # If needed, second attempt to get assembly ID
                if [[ ${#assemblies[@]} -eq 0 ]]; then
                    mapfile -t assemblies < <(esearch -db "$dl_db" -query "$accession" | elink -target assembly |
                        esummary | xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)
                fi
            else
                # 'WP_' (multispecies) entries
                mapfile -t assemblies < <(esearch -db "$dl_db" -query "$accession" |
                    elink -target nuccore -name protein_nuccore_wp |
                    elink -db nuccore -target assembly -name nuccore_assembly |
                    esummary | xtract -pattern DocumentSummary -element AssemblyAccession)
            fi
        fi
        
        # Add the retrieved assemblies to the assembly list & lookup file
        if [[ ${#assemblies[@]} -gt 0 ]]; then
            log_time "Found ${#assemblies[@]} assemblies for accession $accession" 
            echo "${assemblies[@]}" | tr " " "\n" > "$assembly_list_acc"
            cat "$assembly_list_acc" >> "$assembly_list"

            awk -v accession="$accession" '{print $0 "\t" accession}' "$assembly_list_acc" | tee -a "$assembly_lookup"
            #echo -e "${accession}\t${assembly}" | tee -a "$assembly_lookup"
            rm "$assembly_list_acc"
        else
            log_time "WARNING: No assembly found for subject $accession"
        fi
    done 9< <(cut -f 2 "$blast_out_final" | sort -u)
    
    # Report
    log_time "Listing the assembly list and subject-to-assembly lookup table files:"
    ls -lh "$assembly_list" "$assembly_lookup"
    log_time "Number of distinct genomes to be downloaded: $(wc -l < "$assembly_list")"

    # Download genome metadata
    log_time "Getting the genome metadata..."
    runstats datasets summary genome accession \
        --inputfile "$assembly_list" --as-json-lines |
        dataformat tsv genome --fields "$META_FIELDS" > "$meta_file"
    log_time "Listing the metadata file..."
    ls -lh "$meta_file"
}

dl_genomes() {
    log_time "Now downloading full genomes of matched sequences..."

    # Download genomes
    runstats datasets download genome accession \
        --inputfile "$assembly_list" \
        --include genome \
        --filename "$outdir"/genomes/genomes.zip \
        --no-progressbar \
        --api-key "$NCBI_API_KEY"

    # Unzip the genomes ZIP file
    unzip -q -o "$outdir"/genomes/genomes.zip -d "$outdir"/genomes

    # Rename genome files
    while IFS= read -r -d '' file; do
        acc_nr=$(dirname "$file" | xargs basename)
        outfile="$outdir"/genomes/"$acc_nr"."${file##*.}"
        mv "$file" "$outfile"
    done < <(find "$outdir"/genomes/ncbi_dataset -type f -wholename "*data/GC*" -print0)

    # Clean up & report
    rm -r "$outdir"/genomes/README.md "$outdir"/genomes/ncbi_dataset "$outdir"/genomes/genomes.zip
    log_time "Listing the downloaded genomes..."
    ls -lh "$outdir"/genomes/*fna
}

dl_subjects() {
    echo -e "\n================================================================"
    log_time "Now downloading full subjects..."
    log_time "Number of downloads: $(cut -f 2 "$blast_out_final" | sort -u | wc -l)"
    mkdir -p "$outdir"/subjects/concat

    while read -r accession; do
        log_time "Subject: $accession"
        outfile="$outdir"/subjects/"$accession".fa
        efetch -db "$dl_db" -format fasta -id "$accession" > "$outfile"
    done < <(cut -f 2 "$blast_out_final" | sort -u)

    log_time "Listing the subject output FASTA files:"
    ls -lh "$outdir"/subjects

    log_time "Creating a multi-FASTA file with all subject sequences..."
    cat "$outdir"/subjects/*fa > "$outdir"/subjects/concat/all.fa
    ls -lh "$outdir"/subjects/concat/all.fa
}

dl_aligned() {
    echo -e "\n================================================================"
    log_time "Now downloading the aligned parts of subjects..."
    log_time "Number of downloads: $(cut -f 2,9,10 "$blast_out_final" | sort -u | wc -l)"
    mkdir -p "$outdir"/aligned/concat

    while read -r accession start stop; do
        log_time "Subject: $accession     Start pos: $start     Stop pos: $stop"
        outfile="$outdir"/aligned/"$accession"_"$start"-"$stop".fa
        efetch -db "$dl_db" -format fasta \
            -id "$accession" -seq_start "$start" -seq_stop "$stop" > "$outfile"
    done < <(cut -f 2,9,10 "$blast_out_final" | sort -u)

    log_time "Listing the aligned-only output files:"
    ls -lh "$outdir"/aligned

    log_time "Creating a multi-FASTA file with all aligned-only sequences..."
    cat "$outdir"/aligned/*fa > "$outdir"/aligned/concat/all.fa
    ls -lh "$outdir"/aligned/concat/all.fa
}

dl_nuc_from_prot() {
    # NOTE: This won't attempt to download "WP_" accessions, since these are
    #       multispecies proteins not directly linked to a nucleotide sequence

    echo -e "\n================================================================"
    log_time "Now downloading the nucleotide sequences for protein hits..."
    mkdir -p "$outdir"/nuc_from_prot

    while read -r accession; do
        echo "Subject: $accession"
        outfile="$outdir"/nuc_from_prot/"$accession".fa
        efetch -db protein -format fasta_cds_na -id "$accession" > "$outfile"
    done < <(cut -f 2 "$blast_out_final" | sort -u | grep -v "WP_")
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
outdir=
task_opt=
tax_ids= && tax_opt= && tax_optarg=
max_target_seqs= && maxtarget_opt=
spacer=
db=
threads= && thread_opt=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --out_format )      shift && blast_format=$1 ;;
        --force )           force=true ;;
        --no_header )       add_header=false ;;
        --max_target_seqs ) shift && max_target_seqs=$1 ;;
        --tax_ids )         shift && tax_ids=$1 ;;
        --db )              shift && db=$1 ;;
        --blast_type )      shift && blast_type=$1 ;;
        --blast_task )      shift && blast_task=$1 ;;
        --local )           local=true && remote_opt= ;;
        --top_n_query )     shift && top_n_query=$1 ;;
        --top_n_subject )   shift && top_n_subject=$1 ;;
        --evalue )          shift && evalue=$1 ;;
        --pct_id )          shift && pct_id=$1 ;;
        --pct_cov )         shift && pct_cov=$1 ;;
        --find_genomes )    to_find_genomes=true ;;
        --dl_genomes )      to_dl_genomes=true ;;
        --dl_subjects )     to_dl_subjects=true ;;
        --dl_aligned )      to_dl_aligned=true ;;
        --env )             shift && env=$1 ;;
        --no_strict )       strict_bash=false ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        --version )         version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"

# Make paths absolute
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"
[[ ! "$infile" =~ ^/ ]] && infile="$PWD"/"$infile"

# Define output files
blast_out_raw="$outdir"/blast_out_raw.tsv
blast_out_sorted="$outdir"/blast_out_sorted.tsv
blast_out_final="$outdir"/blast_out_final.tsv

# Create the output dirs
mkdir -p "$outdir"/logs

# If download-genomes is true, so should find-genomes be
[[ "$to_dl_genomes" == "true" ]] && to_find_genomes=true

# BLAST task option
[[ -n "$blast_task" ]] && task_opt=" -task $blast_task"

# Build BLAST taxon ID option (format: -entrez_query "txid343[Organism:exp] OR txid56448[Organism:exp]")
if [[ -n "$tax_ids" ]]; then
    spacer=" "
    tax_opt=" -entrez_query"
    IFS=',' read -ra tax_array <<< "$tax_ids"

    for tax_id in "${tax_array[@]}"; do
        if [[ -z "$tax_optarg" ]]; then
            tax_optarg="txid$tax_id[Organism:exp]"
        else
            tax_optarg="$tax_optarg OR txid$tax_id[Organism:exp]"
        fi
    done
fi

# Max target sequences
[[ -n "$max_target_seqs" ]] && maxtarget_opt=" -max_target_seqs $max_target_seqs"

# Defaults based on the type of BLAST
[[ "$blast_type" == "blastx" || "$blast_type" == "blastp" ]] && db_type=prot
[[ "$db_type" == "prot" ]] && dl_db=protein

if [[ "$local" == true ]]; then
    # Local BLAST nr of cores
    set_threads "$IS_SLURM"
    thread_opt=" -num_threads $threads"

    # Local BLAST default db
    if [[ "$db_type" == "prot" ]]; then
        [[ -z "$db" ]] && db="$DEFAULT_LOCAL_DB_AA"
    else
        [[ -z "$db" ]] && db="$DEFAULT_LOCAL_DB_NT"
    fi
else
    # Remote BLAST default DB
    if [[ "$db_type" == "prot" ]]; then
        [[ -z "$db" ]] && db="$remote_db_aa"
    else
        [[ -z "$db" ]] && db="$remote_db_nt"
    fi
fi

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Output dir:                               $outdir"
echo
echo "BLAST db:                                 $db"
echo "BLAST type:                               $blast_type"
[[ -n "$blast_task" ]] && echo "BLAST task:                               $blast_task"
echo "Run BLAST locally?                        $local"
echo
echo "Force BLAST run even if output exists?    $force"
echo "Add column header to BLAST output?        $add_header"
echo
echo "Evalue threshold:                         $evalue"
[[ -n "$pct_id" ]] && echo "Percent identity threshold:               $pct_id"
[[ -n "$pct_cov" ]] && echo "Alignment coverage threshold:             $pct_cov"
[[ -n "$top_n_query" ]] && echo "Filter to top N hits per query:           $top_n_query"
[[ -n "$top_n_subject" ]] && echo "Filter to top N hits per subject:         $top_n_subject"
[[ -n "$max_target_seqs" ]] && echo "Max. nr. of target sequences:             $max_target_seqs"
if [[ -n "$tax_ids" ]]; then
    echo
    echo "Taxon IDs:                                $tax_ids"
    echo "Nr of taxonomic IDs:                      ${#tax_array[@]}"
    echo "Taxonomic selection option:               $tax_optarg"
fi
echo
echo "Find genome accessions?                   $to_find_genomes"
echo "Download full genomes?                    $to_dl_genomes"
echo "Download full subjects?                   $to_dl_subjects"
echo "Download aligned parts of sequences?      $to_dl_aligned"
echo
echo "Number of queries in the input file:      $(grep -c "^>" "$infile")"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                              RUN
# ==============================================================================
# Run BLAST
if [[ -s "$blast_out_raw" && "$force" == false ]]; then
    log_time "Skipping BLAST, output file exists ($blast_out_raw) and --force is false..."
else
    run_blast
fi

# Process BLAST output
process_blast
[[ "$n_subjects" -eq 0 ]] && echo "EXITING: No BLAST hits remained after filtering" && exit 0

# Download aligned parts of sequences
[[ "$to_dl_aligned" == true ]] && dl_aligned

# Download full subjects (matched contigs/Genbank entries only)
[[ "$to_dl_subjects" == true ]] && dl_subjects

# Download full genomes
[[ "$to_find_genomes" == true ]] && find_genomes
[[ "$to_dl_genomes" == true ]] && dl_genomes

# Download gene nucleotide sequences for protein hits
[[ "$to_dl_subjects" == true && "$db_type" == "prot" ]] && dl_nuc_from_prot

# Add header to the final BLAST output file
if [[ "$add_header" == true ]]; then
    header=$(echo "$blast_format" | sed 's/6 //' | tr " " "\t") 
    sed -i "1s/^/$header\n/" "$blast_out_final"
fi

# Final reporting
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"


# ==============================================================================
#                              SANDBOX
# ==============================================================================
#? Include organism taxonomy
## First, create a lookup table with accession numbers and organism taxonomy
## (See https://www.biostars.org/p/367121/, Accession nr examples for testing: accession=KF772785.1 / accession=MH231153.1)
# while read -r line; do
#    accession=$(echo "$line" | cut -f 2)
#    esearch_result=$(esearch -db nuccore -query "$accession" </dev/null)
#    
#    # Get the taxon name
#    taxon=$(echo "$esearch_result" |
#                elink -target taxonomy |
#                efetch -format native -mode xml |
#                grep "ScientificName" | head -n1 |
#                sed -E 's@</?ScientificName>@@g' | sed -e 's/^[ \t]*//')
#
#    echo -e "${accession}\t${taxon}"
#    echo -e "${accession}\t${taxon}" >&2
#    sleep 5s
# done < "$blast_out_top" > "$organism_lookup"
## Second, merge the BLAST output table with the taxonomy lookup table
# join -t $'\t' -1 2 -2 1 "$blast_out_top" "$organism_lookup" > "$blast_out_proc"

#? Alternative way to get aligned sequences - faster
#blastn -remote -task blastn -db nt -outfmt '6 qseqid sacc sseq' -evalue 1e-6 \
#    -query "$query_fa" -out blast_seq.tsv
#awk -v OFS="\n" '{print ">"$2, $3}' blast_seq.tsv > blast_seqs.fa
