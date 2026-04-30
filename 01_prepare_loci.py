#!/usr/bin/env python3
"""
01_prepare_loci.py
==================
Parse the Ulmus microsatellite CSV and emit:
  - loci.bed         : BED4 of the STR region on each reference contig
  - primers.tsv      : per-locus forward / reverse primer sequences
  - loci_config.tsv  : full metadata used by downstream scripts
  - str_regions.bed  : BED for samtools view (amplicon window with padding)

Usage:
    python3 01_prepare_loci.py \
        --csv  Ulmus_amplicon_targets.csv \
        --outdir loci/

Requirements: pandas
"""

import argparse
import os
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv",    required=True,  help="Input CSV file")
    p.add_argument("--outdir", default="loci", help="Output directory")
    p.add_argument("--padding", type=int, default=50,
                   help="bp of padding around STR for samtools fetch (default 50)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)

    # Drop trailing NaN row
    df = df.dropna(subset=["gene"])

    # Rename for convenience
    df = df.rename(columns={
        "gene":                                    "locus",
        "asssay_seq_name":                         "ref_contig",
        "forward_primer_sequence":                 "fwd_primer",
        "reverse_primer(but forward)":             "rev_primer_fwd_orient",
        "reverse_primer_sequence (reverse_complement)": "rev_primer_rc",
        "repeat_motif":                            "motif",
        "motif_length":                            "motif_len",
        "number_repeats":                          "n_repeats",
        "length_of_STR":                           "str_len_ref",
        "Start":                                   "str_start",   # 1-based in CSV
        "Stop":                                    "str_stop",
        "length_assay":                            "amplicon_len",
    })

    # Validate required columns
    required = ["locus","ref_contig","fwd_primer","rev_primer_fwd_orient",
                "motif","motif_len","str_start","str_stop"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # Convert to 0-based half-open BED coordinates
    df["bed_start"] = df["str_start"].astype(int) - 1
    df["bed_stop"]  = df["str_stop"].astype(int)

    # Padded window for read retrieval
    df["fetch_start"] = (df["bed_start"] - args.padding).clip(lower=0)
    df["fetch_stop"]  = df["bed_stop"]  + args.padding

    # ---- loci.bed (STR region only) ----------------------------------------
    bed_path = os.path.join(args.outdir, "loci.bed")
    df[["ref_contig","bed_start","bed_stop","locus"]].to_csv(
        bed_path, sep="\t", header=False, index=False
    )
    print(f"[01] Wrote {bed_path}  ({len(df)} loci)")

    # ---- str_regions.bed (padded, for samtools view) -----------------------
    fetch_bed_path = os.path.join(args.outdir, "str_regions_padded.bed")
    df[["ref_contig","fetch_start","fetch_stop","locus"]].to_csv(
        fetch_bed_path, sep="\t", header=False, index=False
    )
    print(f"[01] Wrote {fetch_bed_path}")

    # ---- primers.tsv -------------------------------------------------------
    primer_path = os.path.join(args.outdir, "primers.tsv")
    df[["locus","fwd_primer","rev_primer_fwd_orient","rev_primer_rc"]].to_csv(
        primer_path, sep="\t", index=False
    )
    print(f"[01] Wrote {primer_path}")

    # ---- loci_config.tsv (full metadata) -----------------------------------
    config_cols = ["locus","ref_contig","fwd_primer","rev_primer_fwd_orient",
                   "rev_primer_rc","motif","motif_len","str_len_ref",
                   "bed_start","bed_stop","fetch_start","fetch_stop","amplicon_len"]
    # motif_len may be float from CSV; cast safely
    df["motif_len"] = df["motif_len"].astype(int)
    df["str_len_ref"] = df["str_len_ref"].fillna(0).astype(int)

    config_path = os.path.join(args.outdir, "loci_config.tsv")
    df[config_cols].to_csv(config_path, sep="\t", index=False)
    print(f"[01] Wrote {config_path}")

    # ---- summary -----------------------------------------------------------
    print(f"\n[01] Summary")
    print(f"     Loci: {len(df)}")
    print(f"     Motif lengths: {sorted(df['motif_len'].unique())}")
    print(f"     Median STR length (ref): {df['str_len_ref'].median()} bp")
    print(f"     Amplicon sizes range: {int(df['amplicon_len'].min())}–{int(df['amplicon_len'].max())} bp")
    print(f"     NOTE: {(df['amplicon_len'] > 150).sum()} loci have amplicons > 150 bp "
          f"(span Illumina reads — handled by alignment + region-level pileup)")


if __name__ == "__main__":
    main()
