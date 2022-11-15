#!/bin/bash

## Bash strict mode
set -euo pipefail

## Parse command-line args
list_in=$1
yaml_out=$2

## Create template file
cat > "$yaml_out".tmp1 <<'_EOF'
[
    {
        orientation: "fr",
        type: "paired-end",
        right reads: [
        R2_reads
        ],
        left reads: [
        R1_reads
        ]
    },
]
_EOF

## Make dirs absolute if needed
sed "/^\//d" "$list_in" | sed -E "s@^@$PWD/@" > "$yaml_out".modlist
sed -n "/^\//p" "$list_in" >> "$yaml_out".modlist   # Paths that were already absolute

## Replace placeholder strings with filenames
R1=$(sed -e 's/^/"/' -e 's/$/"/' "$yaml_out".modlist | sed 's/$/,/')
R2=$(sed -e 's/^/"/' -e 's/$/"/' "$yaml_out".modlist | sed 's/$/,/' | sed 's/_R1/_R2/')
awk -v r="$R1" '{gsub(/R1_reads/,r)}1' "$yaml_out".tmp1 > "$yaml_out".tmp2
awk -v r="$R2" '{gsub(/R2_reads/,r)}1' "$yaml_out".tmp2 > "$yaml_out"

## Remove temporay files
rm "$yaml_out".tmp1 "$yaml_out".tmp2 "$yaml_out".modlist

## Print output file to screen
cat "$yaml_out"