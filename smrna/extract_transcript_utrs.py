#!/usr/bin/env python3
"""Extract transcript-level 3' UTR sequences from a GTF/GFF using bedtools.

- Run with `-h` flag for usage instructions.
- NOTE: `bedtools` must be installed and in the PATH.

Outputs:
- Transcript-level UTR FASTA (usable by Miranda)
- Optional BED12 intermediate (for reproducibility)
- Optional TargetScan-style UTR table: transcript_id<TAB>species_id<TAB>sequence
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


def parse_attrs(attr: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
            continue
        fields = part.split(None, 1)
        if len(fields) == 2:
            out[fields[0].strip()] = fields[1].strip().strip('"')
    return out


def get_tx_id(attrs: Dict[str, str]) -> str:
    tx = attrs.get("transcript_id", "")
    if tx:
        return tx
    parent = attrs.get("Parent", "")
    if parent.startswith("transcript:"):
        parent = parent[len("transcript:") :]
    return parent


def run(cmd: List[str]) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"ERROR: command failed ({exc.returncode}): {' '.join(cmd)}")


def normalize_fasta_headers(path: str) -> None:
    tmp = path + ".tmp"
    with (
        open(path, "r", encoding="utf-8") as fin,
        open(tmp, "w", encoding="utf-8") as fout,
    ):
        for line in fin:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                name = re.sub(r"\([+-]\)$", "", name)
                fout.write(f">{name}\n")
            else:
                fout.write(line)
    os.replace(tmp, path)


def write_targetscan_table(fasta_path: str, out_path: str, species_id: str) -> None:
    with (
        open(fasta_path, "r", encoding="utf-8") as fin,
        open(out_path, "w", encoding="utf-8") as fout,
    ):
        tx = None
        seq_parts: List[str] = []
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if tx is not None:
                    fout.write(f"{tx}\t{species_id}\t{''.join(seq_parts)}\n")
                tx = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if tx is not None:
            fout.write(f"{tx}\t{species_id}\t{''.join(seq_parts)}\n")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Extract transcript-level 3' UTR FASTA from GTF/GFF + genome."
    )
    p.add_argument("--gtf", required=True, help="Annotation file (GTF or GFF)")
    p.add_argument("--genome", required=True, help="Genome FASTA")
    p.add_argument(
        "--out-fasta", required=True, help="Output transcript-level UTR FASTA"
    )
    p.add_argument("--out-bed12", default=None, help="Optional output BED12 path")
    p.add_argument(
        "--out-targetscan", default=None, help="Optional output TargetScan UTR table"
    )
    p.add_argument(
        "--species-id", default="9999", help="Species ID for --out-targetscan"
    )
    p.add_argument(
        "--utr-feature",
        default="three_prime_utr",
        help="Feature name to extract from column 3",
    )
    args = p.parse_args()

    tx_parts: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)

    with open(args.gtf, "r", encoding="utf-8") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feature, start_s, end_s, _, strand, _, attr_s = cols
            if feature != args.utr_feature:
                continue
            attrs = parse_attrs(attr_s)
            tx = get_tx_id(attrs)
            if not tx:
                continue
            start0 = int(start_s) - 1
            end = int(end_s)
            tx_parts[tx].append((chrom, start0, end, strand))

    if not tx_parts:
        raise SystemExit(
            f"ERROR: no '{args.utr_feature}' features with transcript IDs found in {args.gtf}"
        )

    bed12_path = args.out_bed12 or (args.out_fasta + ".bed12")

    with open(bed12_path, "w", encoding="utf-8") as out_bed:
        for tx in sorted(tx_parts):
            parts = tx_parts[tx]
            chroms = {c for c, _, _, _ in parts}
            strands = {s for _, _, _, s in parts}
            if len(chroms) != 1 or len(strands) != 1:
                continue

            chrom = next(iter(chroms))
            strand = next(iter(strands))
            parts_sorted = sorted(parts, key=lambda x: x[1])

            chrom_start = min(s for _, s, _, _ in parts_sorted)
            chrom_end = max(e for _, _, e, _ in parts_sorted)

            block_sizes = [str(e - s) for _, s, e, _ in parts_sorted]
            block_starts = [str(s - chrom_start) for _, s, _, _ in parts_sorted]

            row = [
                chrom,
                str(chrom_start),
                str(chrom_end),
                tx,
                "0",
                strand,
                str(chrom_start),
                str(chrom_end),
                "0",
                str(len(parts_sorted)),
                ",".join(block_sizes) + ",",
                ",".join(block_starts) + ",",
            ]
            out_bed.write("\t".join(row) + "\n")

    run(
        [
            "bedtools",
            "getfasta",
            "-fi",
            args.genome,
            "-bed",
            bed12_path,
            "-split",
            "-s",
            "-nameOnly",
            "-fo",
            args.out_fasta,
        ]
    )

    normalize_fasta_headers(args.out_fasta)

    if args.out_targetscan:
        write_targetscan_table(args.out_fasta, args.out_targetscan, args.species_id)


if __name__ == "__main__":
    main()
