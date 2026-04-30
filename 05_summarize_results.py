#!/usr/bin/env python3
"""
05_summarize_results.py
=======================
Merge all calls, apply final QC thresholds, and write the canonical output:

  results/final_alleles.csv      — primary deliverable: one row per sample-locus
  results/genotype_matrix.csv    — wide pivot for downstream popgen tools
  results/locus_stats.csv        — per-locus summary (call rate, allelic diversity)
  results/sample_stats.csv       — per-sample summary (call rate, mean coverage)
  results/pipeline_report.txt    — human-readable summary

The genotype_matrix.csv uses the format compatible with:
  - polysat (R)         : alleles as "120/120/132/132" (tetraploid)
  - pegas (R)           : same slash notation
  - Structure/STRUCTURE : will need manual conversion from the alleles columns

Usage:
    python3 05_summarize_results.py \
        --calls    calls/allele_calls.tsv \
        --config   loci/loci_config.tsv \
        --outdir   results/ \
        [--min-call-rate 0.7] [--min-locus-call-rate 0.5]

Requirements: pandas, numpy
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls",                required=True)
    p.add_argument("--config",               required=True)
    p.add_argument("--outdir",               default="results")
    p.add_argument("--min-call-rate",        type=float, default=0.70,
                   help="Min fraction of loci called per sample (default 0.70)")
    p.add_argument("--min-locus-call-rate",  type=float, default=0.50,
                   help="Min fraction of samples called per locus (default 0.50)")
    p.add_argument("--missing-code",         default="./.",
                   help="Code used for missing genotypes (default ./. )")
    return p.parse_args()


def allele_richness(genotypes, missing_code):
    """Count distinct allele lengths across all genotypes for a locus."""
    all_alleles = set()
    for g in genotypes:
        if g == missing_code or pd.isna(g):
            continue
        parts = str(g).split("/")
        for a in parts:
            if a not in (".", ""):
                all_alleles.add(a)
    return len(all_alleles)


def n_alleles_per_sample(genotype, missing_code):
    if genotype == missing_code or pd.isna(genotype):
        return np.nan
    alleles = [a for a in str(genotype).split("/") if a not in (".", "")]
    return len(set(alleles))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    calls  = pd.read_csv(args.calls, sep="\t")
    config = pd.read_csv(args.config, sep="\t")

    print(f"[05] Loaded {len(calls)} call records")
    print(f"[05] Samples: {calls['sample'].nunique()}")
    print(f"[05] Loci:    {calls['locus'].nunique()}")

    mc = args.missing_code
    calls["is_called"] = calls["genotype"] != mc

    # ---- per-locus stats --------------------------------------------------
    locus_stats_rows = []
    for locus, grp in calls.groupby("locus"):
        n_total   = len(grp)
        n_called  = grp["is_called"].sum()
        call_rate = n_called / n_total if n_total > 0 else 0
        ar        = allele_richness(grp["genotype"], mc)
        mean_cov  = grp["total_reads"].mean()

        # Motif info from config
        cfg = config[config["locus"] == locus]
        motif     = cfg["motif"].iloc[0] if len(cfg) > 0 else "?"
        motif_len = cfg["motif_len"].iloc[0] if len(cfg) > 0 else "?"

        locus_stats_rows.append({
            "locus":          locus,
            "motif":          motif,
            "motif_len":      motif_len,
            "n_samples":      n_total,
            "n_called":       n_called,
            "call_rate":      round(call_rate, 3),
            "allele_richness": ar,
            "mean_coverage":  round(mean_cov, 1),
            "pass_filter":    call_rate >= args.min_locus_call_rate,
        })

    locus_stats = pd.DataFrame(locus_stats_rows)
    passing_loci = set(locus_stats[locus_stats["pass_filter"]]["locus"])
    locus_stats_path = os.path.join(args.outdir, "locus_stats.csv")
    locus_stats.to_csv(locus_stats_path, index=False)
    print(f"[05] Wrote {locus_stats_path}")
    print(f"[05] Loci passing call-rate filter: {len(passing_loci)}/{len(locus_stats)}")

    # ---- per-sample stats -------------------------------------------------
    sample_stats_rows = []
    for sample, grp in calls.groupby("sample"):
        n_total   = len(grp)
        n_called  = grp["is_called"].sum()
        call_rate = n_called / n_total if n_total > 0 else 0
        ploidy    = grp["ploidy"].iloc[0]
        mean_cov  = grp["total_reads"].mean()

        sample_stats_rows.append({
            "sample":         sample,
            "ploidy":         ploidy,
            "n_loci_total":   n_total,
            "n_loci_called":  n_called,
            "call_rate":      round(call_rate, 3),
            "mean_coverage":  round(mean_cov, 1),
            "pass_filter":    call_rate >= args.min_call_rate,
        })

    sample_stats = pd.DataFrame(sample_stats_rows)
    passing_samples = set(sample_stats[sample_stats["pass_filter"]]["sample"])
    sample_stats_path = os.path.join(args.outdir, "sample_stats.csv")
    sample_stats.to_csv(sample_stats_path, index=False)
    print(f"[05] Wrote {sample_stats_path}")
    print(f"[05] Samples passing call-rate filter: {len(passing_samples)}/{len(sample_stats)}")

    # ---- final_alleles.csv (primary deliverable) --------------------------
    # Retain only passing samples and passing loci
    final = calls[
        calls["sample"].isin(passing_samples) &
        calls["locus"].isin(passing_loci)
    ].copy()

    # Expand alleles into individual columns for clarity
    # Determine max ploidy in dataset
    max_ploidy = int(final["ploidy"].max()) if len(final) > 0 else 2

    def split_allele(genotype, idx):
        if genotype == mc or pd.isna(genotype):
            return None
        parts = [a for a in str(genotype).split("/") if a not in (".", "")]
        return int(parts[idx]) if idx < len(parts) else None

    for i in range(max_ploidy):
        final[f"allele_{i+1}"] = final["genotype"].apply(lambda g: split_allele(g, i))

    final_cols = (
        ["sample","locus","ploidy","genotype"] +
        [f"allele_{i+1}" for i in range(max_ploidy)] +
        ["dosages","n_peaks","total_reads","peak_fractions","flag"]
    )
    final_cols = [c for c in final_cols if c in final.columns]

    final_path = os.path.join(args.outdir, "final_alleles.csv")
    final[final_cols].sort_values(["sample","locus"]).to_csv(final_path, index=False)
    print(f"[05] Wrote {final_path}  "
          f"({final['sample'].nunique()} samples × {final['locus'].nunique()} loci)")

    # ---- genotype_matrix.csv (wide pivot) ---------------------------------
    geno_wide = final.pivot_table(
        index="sample", columns="locus", values="genotype", aggfunc="first"
    )
    geno_wide.columns.name = None
    geno_wide = geno_wide.fillna(mc)
    geno_path = os.path.join(args.outdir, "genotype_matrix.csv")
    geno_wide.to_csv(geno_path)
    print(f"[05] Wrote {geno_path}")

    # ---- pipeline_report.txt ---------------------------------------------
    report_path = os.path.join(args.outdir, "pipeline_report.txt")
    with open(report_path, "w") as fh:
        fh.write("=" * 60 + "\n")
        fh.write("Ulmus STR Pipeline — Summary Report\n")
        fh.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write("=" * 60 + "\n\n")

        fh.write("INPUT\n")
        fh.write(f"  Total samples:       {calls['sample'].nunique()}\n")
        fh.write(f"  Total loci:          {calls['locus'].nunique()}\n")
        fh.write(f"  Ploidy breakdown:\n")
        for ploidy, cnt in calls.groupby("ploidy")["sample"].nunique().items():
            fh.write(f"    {ploidy}x: {cnt} samples\n")
        fh.write("\n")

        fh.write("AFTER QC FILTERS\n")
        fh.write(f"  Min sample call rate: {args.min_call_rate}\n")
        fh.write(f"  Min locus call rate:  {args.min_locus_call_rate}\n")
        fh.write(f"  Samples passing:      {len(passing_samples)}/{calls['sample'].nunique()}\n")
        fh.write(f"  Loci passing:         {len(passing_loci)}/{calls['locus'].nunique()}\n\n")

        fh.write("LOCUS DETAILS\n")
        fh.write(locus_stats.sort_values("call_rate", ascending=False).to_string(index=False))
        fh.write("\n\n")

        fh.write("SAMPLE DETAILS\n")
        fh.write(sample_stats.sort_values("call_rate", ascending=False).to_string(index=False))
        fh.write("\n")

    print(f"[05] Wrote {report_path}")
    print(f"\n[05] ✓ Pipeline complete. Primary output: {final_path}")
    print(f"     Genotype matrix for popgen: {geno_path}")


if __name__ == "__main__":
    main()
