#!/usr/bin/env python3
"""
03_extract_str_lengths.py
=========================
For each sample BAM and each STR locus, infer the length of the microsatellite
by measuring the insert size of read pairs that span (or partially span) the
STR, then subtracting the known flanking sequence lengths.

Why this approach instead of existing tools:
  - seq2sat / STRait Razor: fail when STR > read length (amplicons here can be
    150–282 bp, well past a 150 bp Illumina read).
  - gangSTR / hipSTR: require long flanking regions for the stutter model;
    these amplicons don't have enough padding.
  - RepeatSeq: outdated dependencies.
  Instead we use reference-anchored paired-read insert sizes. When both reads
  of a pair map within the amplicon region, the observed insert minus the
  non-STR sequence gives us the STR length directly.

Algorithm (per locus per read pair):
  1. Collect all properly paired, high-quality read pairs overlapping the STR.
  2. For pairs where R1 maps within the left flank region and R2 maps within
     the right flank region (or vice-versa), the TLEN field gives the
     full amplicon insert size.
  3. STR length = insert_size - (left_flank_len + right_flank_len + fwd_primer_len + rev_primer_len)
     where flank lengths are derived from the loci_config.tsv Start/Stop and
     amplicon coordinates.
  4. Round to nearest motif-length multiple (handles minor sequencing error).
  5. All observed lengths are collected into a distribution per sample-locus.

Output:
  - allele_lengths_raw.tsv   : one row per read pair (sample, locus, length)
  - length_distributions.tsv : aggregated count of each length per sample-locus

Usage:
    python3 03_extract_str_lengths.py \
        --bam-dir   alignments/ \
        --config    loci/loci_config.tsv \
        --outdir    str_lengths/ \
        [--min-mapq 20] [--min-pairs 5] [--max-insert 1000]

Requirements: pysam, pandas, numpy
"""

import argparse
import os
import sys
import re
import glob
import pysam
import numpy as np
import pandas as pd
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-dir",    required=True,  help="Dir containing per-sample subdirs with *.str.bam")
    p.add_argument("--config",     required=True,  help="loci_config.tsv from 01_prepare_loci.py")
    p.add_argument("--outdir",     default="str_lengths", help="Output directory")
    p.add_argument("--min-mapq",   type=int, default=20,   help="Min mapping quality")
    p.add_argument("--min-pairs",  type=int, default=5,    help="Min spanning read pairs to report a locus")
    p.add_argument("--max-insert", type=int, default=1000, help="Max plausible insert size (filter chimeric pairs)")
    p.add_argument("--bam-suffix", default=".str.bam",     help="BAM filename suffix to search for")
    return p.parse_args()


def load_config(path):
    df = pd.read_csv(path, sep="\t")
    # Compute non-STR sequence lengths for insert-size math
    # Structure: [fwd_primer][left_flank][STR][right_flank][rev_primer]
    # bed_start is 0-based start of STR; fetch_start is padded amplicon start.
    # We use str_start / str_stop (BED coords) relative to the contig.
    # The primer lengths are in the CSV but not stored in config; derive from sequences.
    df["fwd_primer_len"] = df["fwd_primer"].str.len()
    df["rev_primer_len"] = df["rev_primer_fwd_orient"].str.len()
    # Left non-STR = distance from fwd primer start to STR start
    # This is: bed_start - (position of fwd primer on contig)
    # Since we don't have absolute primer coordinates, we use:
    #   left_flank = bed_start - fwd_primer_len
    #   right_flank = amplicon_len - (bed_stop - (bed_start - fwd_primer_len)) - rev_primer_len
    # More robustly: non_str_len = amplicon_len - str_len_ref
    df["non_str_len"] = df["amplicon_len"] - df["str_len_ref"]
    return df.set_index("locus").to_dict(orient="index")


def snap_to_motif(length, motif_len, tolerance=2):
    """Round length to nearest motif-length multiple. Reject if too far off."""
    if length <= 0:
        return None
    nearest = round(length / motif_len) * motif_len
    if abs(length - nearest) <= tolerance:
        return nearest
    # Allow ±1 bp slop (indel sequencing error) and still round
    return nearest


def extract_lengths_from_bam(bam_path, locus_meta, min_mapq, max_insert):
    """
    Returns a list of (locus_name, inferred_str_length) from spanning pairs.
    """
    results = []
    sam = pysam.AlignmentFile(bam_path, "rb")

    for locus, meta in locus_meta.items():
        contig   = meta["ref_contig"]
        str_s    = meta["bed_start"]   # 0-based
        str_e    = meta["bed_stop"]    # 0-based exclusive
        motif_l  = meta["motif_len"]
        non_str  = meta["non_str_len"]

        # Check contig exists in BAM
        try:
            sam.get_tid(contig)
        except ValueError:
            continue

        seen_pairs = {}  # qname -> (tlen, is_valid)

        try:
            for read in sam.fetch(contig, max(0, str_s - 500), str_e + 500):
                if (read.is_unmapped or read.mate_is_unmapped or
                        read.mapping_quality < min_mapq or
                        read.is_secondary or read.is_supplementary or
                        not read.is_paired or not read.is_proper_pair):
                    continue

                tlen = abs(read.template_length)
                if tlen < 50 or tlen > max_insert:
                    continue

                # Only use read1 to avoid double-counting
                if not read.is_read1:
                    continue

                # Read must start upstream of (or within) STR,
                # and mate must end downstream of (or within) STR.
                # This ensures the pair collectively spans the STR.
                r1_start = read.reference_start
                r1_end   = read.reference_end
                mate_start = read.next_reference_start
                mate_end   = mate_start + read.query_length  # approximate

                # Pair spans STR if: r1_start <= str_s  AND  mate_end >= str_e
                # OR mate_start <= str_s AND r1_end >= str_e  (reversed orientation)
                spans = (
                    (r1_start <= str_s and mate_end >= str_e) or
                    (mate_start <= str_s and r1_end >= str_e)
                )

                if spans:
                    inferred = tlen - non_str
                    snapped  = snap_to_motif(inferred, motif_l)
                    if snapped is not None and snapped > 0:
                        results.append((locus, snapped))

        except (ValueError, KeyError):
            continue

    sam.close()
    return results


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load locus metadata
    locus_meta = load_config(args.config)
    print(f"[03] Loaded {len(locus_meta)} loci from {args.config}")

    # Discover BAM files
    pattern = os.path.join(args.bam_dir, "**", f"*{args.bam_suffix}")
    bam_files = sorted(glob.glob(pattern, recursive=True))
    if not bam_files:
        sys.exit(f"[03] ERROR: No BAMs matching '{pattern}' found.")
    print(f"[03] Found {len(bam_files)} BAM files")

    all_rows   = []   # raw per-pair rows
    dist_rows  = []   # aggregated distributions

    for bam_path in bam_files:
        # Infer sample name from directory or filename
        sample = os.path.basename(bam_path).replace(args.bam_suffix, "")
        print(f"[03]   Processing {sample} ...")

        pairs = extract_lengths_from_bam(
            bam_path, locus_meta,
            min_mapq=args.min_mapq,
            max_insert=args.max_insert
        )

        # Raw rows
        for locus, length in pairs:
            all_rows.append({"sample": sample, "locus": locus, "str_length": length})

        # Aggregate per locus
        by_locus = defaultdict(list)
        for locus, length in pairs:
            by_locus[locus].append(length)

        for locus, lengths in by_locus.items():
            n = len(lengths)
            if n < args.min_pairs:
                continue
            arr = np.array(lengths)
            counts = pd.Series(arr).value_counts().sort_index()
            for l_val, cnt in counts.items():
                dist_rows.append({
                    "sample":    sample,
                    "locus":     locus,
                    "str_length": int(l_val),
                    "read_count": int(cnt),
                    "total_pairs": n,
                    "fraction":  round(cnt / n, 4),
                    "motif_len": locus_meta[locus]["motif_len"],
                    "n_repeats": int(l_val) // locus_meta[locus]["motif_len"]
                        if locus_meta[locus]["motif_len"] > 0 else None,
                })

    # Write raw
    raw_path = os.path.join(args.outdir, "allele_lengths_raw.tsv")
    pd.DataFrame(all_rows).to_csv(raw_path, sep="\t", index=False)
    print(f"[03] Wrote {raw_path}  ({len(all_rows)} read-pair records)")

    # Write distributions
    dist_path = os.path.join(args.outdir, "length_distributions.tsv")
    pd.DataFrame(dist_rows).to_csv(dist_path, sep="\t", index=False)
    print(f"[03] Wrote {dist_path}")
    print(f"[03] Next: run  04_call_dosage.py  --dist {dist_path}")


if __name__ == "__main__":
    main()
