The script `cutadapt.sh` will run Cutadapt on all FASTQ files in a directory,
removing primers specified by the user.

## Command-line options

- `-i`  Input dir (REQUIRED)
- `-o`  Output dir (REQUIRED)
- `-f`  Forward primer
- `-r`  Reverse primer
- `-p`  Primer file (one line per primer pair, forward and reverse whitespace-separated)
- `-d`  Don't discard sequences with no primers (default: discard)
- `-h`  Print help.

## Primer sequences

You can provide primer sequences to remove in two ways:

1. For a single primer pair, using the command-line options `-f <forward-primer> -r <reverse-primer>`.

2. For any number of primers pairs, using a file with primer sequences specified
   using `-p <primer-file>`.

   The file should contain one primer pair per line, with the forward primer first,
   then a space or tab, and then the reverse primer.

   For example:
   
   ```
   primer1_f primer1_r
   primer2_f primer2_r
   ```
   
   Or with actual sequences:
   ```
   GAGTGYCAGCMGCCGCGGTAA ACGGACTACNVGGGTWTCTAAT
   CAGTGYCAGCMGCCGCGGTAA TCGGACTACNVGGGTWTCTAAT
   ```

Note that you don't need to (and should not) provide the reverse complements
of the primers, since these will be computed by the script.

## An example with a primer pair provided on the command line:

```sh
indir=data/fastq/raw
outdir=data/fastq/trimmed
primer_f=GAGTGYCAGCMGCCGCGGTAA
primer_r=ACGGACTACNVGGGTWTCTAAT
sbatch cutadapt.sh -i "$indir" -o "$outdir" -f "$primer_f" -r "$primer_r"
```

## An example with primers provided in a file:

```sh
indir=data/fastq/raw
outdir=data/fastq/trimmed
primer_file=data/meta/primers.txt
sbatch cutadapt.sh -i "$indir" -o "$outdir" -p "$primer_file"
```