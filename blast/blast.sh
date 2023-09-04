#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=blast
#SBATCH --output=slurm-blast-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run NCBI BLAST on an input (query) FASTA file,
and optionally download aligned sequences and/or genomes
The FASTA file can contain multiple/many sequences,
though it will be quicker to split a multiFASTA file,
and submit a separate job for each single-sequence FASTA file.
Additionally, download sequences (e.g. with --download_genomes)
will not be separated by query sequence"
SCRIPT_VERSION="2023-09-05"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
VERSION_COMMAND="blastn -version; datasets --version"
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
export LC_ALL=C                     # Locale for sorting

# Defaults - generic
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/blast
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                  # When true, just print tool & script version info and exit

# Constants - settings
META_FIELDS="accession,assminfo-name,organism-name,assminfo-refseq-category,assminfo-level,assmstats-number-of-contigs,assmstats-contig-n50"
BLAST_FORMAT="6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp stitle staxids"
#> Can't get scientific name of subject seq to be included ('ssciname' / 'sscinames')
#> In addition to the 'qcovhsp' included above, there is also 'qcovs', which will contain the total coverage across all HSPs

# Defaults - settings
local=false && remote_opt="-remote" # Run BLAST locally (own db) or remotely (NCBI's db over the internet)
db=nt                               # BLAST db
blast_type=blastn                   # BLAST type
top_n=                              # Keep the top-N hits only (empty => keep all)
evalue="1e-6"                       # E-value threshold
pct_id=                             # % identity threshold (empty => no threshold)
pct_cov=                            # Threshold for % of query covered by the alignment length (empty => no threshold)
force=false                         # Don't rerun BLAST if the output file already exists
to_dl_genomes=false                 # Download full genomes of subjects?
to_dl_subjects=false                # Download full subjects?
to_dl_aligned=false                 # Download aligned sequences?

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
    echo "  - Limit online BLAST database to specific taxa (using NCBI taxon IDs):"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --tax_ids '343,56448'"
    echo "  - Download aligned parts of sequences, full accessions, and full genomes:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --dl_aligned --dl_subjects --dl_genomes"
    echo "  - Use % identity and query coverage thresholds:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --pct_id 90 --pct_cov 90"
    echo "  - Keep only the best 10 hits per query:"
    echo "      sbatch $0 -i my_seq.fa -o results/blast --top_n 10"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input FASTA file (can contain one or more sequences)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "BLAST OPTIONS:"
    echo "  --local                     Run BLAST with a local database                         [default: $local]"
    echo "  --db                <str>   If running remotely: an NCBI database name like 'nt' or 'nr' [default: $db]"
    echo "                              If running locally: the prefix (dir + db name, no file extensions) of a local BLAST db"
    echo "  --blast_type        <str>   Blast type, e.g. 'blastn' or 'blastp'                   [default: $blast_type]"
    echo "  --force                     Run BLAST even if the output file already exists        [default: $force]"
    echo
    echo "BLAST THRESHOLD AND FILTERING OPTIONS:"
    echo "  --tax_ids           <str>   Comma-separated list of NCBI taxon IDs (just the numbers, no 'txid' prefix)"
    echo "                              The BLAST search will be limited to these taxa          [default: use full database]"
    echo "  --evalue            <num>   E-value threshold in scientific notation                [default: $evalue]"
    echo "                                This option will be passed to BLAST, so even the raw"
    echo "                                BLAST output will not contain hits that do not pass this"
    echo "  --pct_id            <num>   Percentage identity threshold                           [default: none]"
    echo "                                This threshold will be applied after running blast"
    echo "  --pct_cov           <num>   Threshold for % of query covered by the alignment       [default: none]"
    echo "                                This threshold will be applied after running blast"
    echo "  --top_n             <int>   Only keep the top N hits for each query                 [default: keep all]"
    echo "                                This threshold will be applied after running blast"
    echo
    echo "SEQUENCE DOWNLOAD OPTIONS:"
    echo "  --dl_aligned        <str>   Download aligned parts of subject (db) sequences        [default: $to_dl_aligned]"
    echo "  --dl_subjects     <str>   Download full subject (db) sequences that were aligned  [default: $to_dl_subjects]"
    echo "  --dl_genomes        <str>   Download full genomes of sequences that were aligned    [default: $to_dl_genomes]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you'll have to provide one in order to run the script with a container.)"
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
    runstats blastn \
        "${remote_opt}""${thread_opt}" \
        -task "$blast_type" \
        -db "$db" \
        -query "$infile" \
        -out "$blast_out_raw" \
        -evalue $evalue \
        -outfmt "$BLAST_FORMAT" \
        "${tax_arg[@]}"
}

process_blast() {
    log_time "Nr of hits in the raw BLAST output file: $(wc -l < "$blast_out_raw")"

    # Sort by: query name (1), then e-value (11), then bitscore (12), then % identical (3)
    log_time "Sorting BLAST output by goodness of the match"
    sort -k1,1 -k11,11g -k12,12gr -k3,3gr "$blast_out_raw" > "$blast_out_sorted"

    # Filter by percent identity
    if [[ -n "$pct_id" ]]; then
        log_time "Filtering output using a percent identity threshold of $pct_id"
        blast_out_id="$outdir"/blast_out_pctid.tsv
        
        awk -F"\t" -v OFS="\t" -v pct_id="$pct_id" \
            '$3 >= pct_id' "$blast_out_sorted" > "$blast_out_id"
        log_time "Nr of retained hits: $(wc -l < "$blast_out_id")"
    else
        blast_out_id="$blast_out_sorted"
    fi

    # Filter by percent coverage
    if [[ -n "$pct_cov" ]]; then
        log_time "Filtering output using a percent coverage threshold of $pct_cov"
        blast_out_cov="$outdir"/blast_out_cov.tsv

        awk -F"\t" -v OFS="\t" -v pct_cov="$pct_cov" \
            '$14 >= pct_cov' "$blast_out_sorted" > "$blast_out_cov"
        log_time "Nr of retained hits: $(wc -l < "$blast_out_cov")"
    else
        blast_out_cov="$blast_out_id"
    fi

    # Only retain top N matches
    if [[ -n "$top_n" ]]; then
        log_time "Getting the top $top_n hits for each query"
        while read -r query; do
            grep -w -m "$top_n" "$query" "$blast_out_sorted"
        done < <(cut -f1 "$blast_out_sorted" | sort -u) > "$blast_out_final"
    else
        cp "$blast_out_cov" "$blast_out_final"
    fi

    # Report & clean
    [[ -f "$blast_out_cov" ]] && rm "$blast_out_cov"
    [[ -f "$blast_out_id" ]] && rm "$blast_out_id"
    [[ -f "$blast_out_sorted" ]] && rm "$blast_out_sorted"

    log_time "Listing the final BLAST output file:"
    ls -lh "$blast_out_final"

    log_time "Showing the first few lines of the final BLAST output file:"
    head -n 5 "$blast_out_final"

    n_accessions=$(cut -f 2 "$blast_out_final" | sort -u | wc -l)
    log_time "Nr of distinct accessions in the final BLAST output file: $n_accessions"
}

dl_genomes() {
    echo -e "\n================================================================"
    log_time "Now downloading full genomes of matched sequences..."
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
        # First attempt to get assembly ID
        #(Note: 'head' at end because in some cases, multiple assembly versions are associated with an accession)
        assembly=$(esearch -db nuccore -query "$accession" | elink -target assembly |
                   esummary | xtract -pattern DocumentSummary -element RefSeq | head -n 1)
        
        # If needed, second attempt to get assembly ID
        if [[ -z "$assembly" ]]; then
            assembly=$(esearch -db nuccore -query "$accession" | elink -target assembly |
                       esummary | xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)
        fi
        
        if [[ -n "$assembly" ]]; then
            echo "$assembly" >> "$assembly_list"
        else
            log_time "WARNING: No assembly found for accession $accession"
        fi
        echo -e "${accession}\t${assembly}" | tee -a "$assembly_lookup"
    
    done 9< <(cut -f 2 "$blast_out_final" | sort -u)
    
    # Report
    log_time "Listing the assembly list and accession-to-assembly lookup table files:"
    ls -lh "$assembly_list" "$assembly_lookup"
    log_time "Number of distinct genomes to be downloaded: $(wc -l < "$assembly_list")"

    # Download genome metadata
    log_time "Getting the genome metadata..."
    runstats datasets summary genome accession \
        --inputfile "$assembly_list" --as-json-lines |
        dataformat tsv genome --fields "$META_FIELDS" > "$meta_file"
    log_time "Listing the metadata file..."
    ls -lh "$meta_file"

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
    ls -lh "$outdir"/genomes
}

dl_subjects() {
    echo -e "\n================================================================"
    log_time "Now downloading full aligned subjects..."
    log_time "Number of downloads: $(cut -f 2 "$blast_out_final" | sort -u | wc -l)"
    mkdir -p "$outdir"/accessions

    while read -r accession; do
        log_time "Accession: $accession"
        outfile="$outdir"/accessions/"$accession".fa
        efetch -db nuccore -format fasta -id "$accession" > "$outfile"
    done < <(cut -f 2 "$blast_out_final" | sort -u)

    log_time "Listing the subject output FASTA files:"
    ls -lh "$outdir"/accessions
}

dl_aligned() {
    echo -e "\n================================================================"
    log_time "Now downloading the aligned parts of subjects..."
    log_time "Number of downloads: $(cut -f 2,9,10 "$blast_out_final" | sort -u | wc -l)"
    mkdir -p "$outdir"/aligned_only/concat

    while read -r accession start stop; do
        log_time "Accession: $accession     Start pos: $start     Stop pos: $stop"
        outfile="$outdir"/aligned_only/"$accession"_"$start"-"$stop".fa
        efetch -db nuccore -format fasta \
            -id "$accession" -seq_start "$start" -seq_stop "$stop" > "$outfile"
    done < <(cut -f 2,9,10 "$blast_out_final" | sort -u)

    log_time "Listing the aligned-only output files:"
    ls -lh "$outdir"/aligned_only

    log_time "Creating a multi-FASTA file with all aligned-only sequences..."
    cat "$outdir"/aligned_only/*fa > "$outdir"/aligned_only/concat/all.fa
    ls -lh "$outdir"/aligned_only/concat/all.fa
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
outdir=
tax_ids=
threads= && thread_opt=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --tax_ids )         shift && tax_ids=$1 ;;
        --db )              shift && db=$1 ;;
        --blast_type )      shift && blast_type=$1 ;;
        --local )           local=true && remote_opt= ;;
        --force )           force=true ;;
        --top_n )           shift && top_n=$1 ;;
        --evalue )          shift && evalue=$1 ;;
        --pct_id )          shift && pct_id=$1 ;;
        --pct_cov )         shift && pct_cov=$1 ;;
        --dl_genomes )      to_dl_genomes=true ;;
        --dl_subjects )   to_dl_subjects=true ;;
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

# Build BLAST taxon ID option (format: -entrez_query "txid343[Organism:exp] OR txid56448[Organism:exp]")
tax_arg=
if [[ -n "$tax_ids" ]]; then
    IFS=',' read -ra tax_array <<< "$tax_ids"
    for tax_id in "${tax_array[@]}"; do
        [[ -n "$tax_arg" ]] && tax_arg="$tax_arg OR txid$tax_id[Organism:exp]"
        [[ -z "$tax_arg" ]] && tax_arg="txid$tax_id[Organism:exp]"
    done
    tax_arg=(-entrez_query "$tax_arg")
fi

# Build BLAST nr threads option (only allowed for local BLAST)
set_threads "$IS_SLURM"
[[ "$local" == true ]] && thread_opt=" -num_threads $threads"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Output dir:                               $outdir"
echo "BLAST db:                                 $db"
echo "BLAST type:                               $blast_type"
echo "Run BLAST locally?                        $local"
echo "Force BLAST run even if output exists?    $force"
echo
echo "Evalue threshold:                         $evalue"
[[ -n "$pct_id" ]] && echo "Percent identity threshold:               $pct_id"
[[ -n "$pct_cov" ]] && echo "Alignment coverage threshold:             $pct_cov"
[[ -n "$top_n" ]] && echo "Subset output to top N hits:              $top_n"
if [[ -n "$tax_ids" ]]; then
    echo
    echo "Taxon IDs:                                $tax_ids"
    echo "Nr of taxonomic IDs:                      ${#tax_array[@]}"
    echo "Taxonomic selection option:               ${tax_arg[*]}"
fi
echo
echo "Download full genomes?                    $to_dl_genomes"
echo "Download full subjects?                   $to_dl_subjects"
echo "Download aligned parts of sequences?      $to_dl_aligned"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                              RUN
# ==============================================================================
# Run BLAST
if [[ -f "$blast_out_raw" && "$force" == false ]]; then
    log_time "Skipping BLAST, output file exists ($blast_out_raw) and --force is false..."
else
    run_blast
fi

# Process BLAST output
process_blast

# Download aligned parts of sequences
[[ "$to_dl_aligned" == true ]] && dl_aligned

# Download full subjects (matched contigs/Genbank entries only)
[[ "$to_dl_subjects" == true ]] && dl_subjects

# Download full genomes
[[ "$to_dl_genomes" == true ]] && dl_genomes

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
