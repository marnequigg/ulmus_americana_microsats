#!/usr/bin/env python3
"""
04_call_dosage.py
=================
Call microsatellite alleles and dosage for diploid (2x) and tetraploid (4x)
samples from the length distributions produced by 03_extract_str_lengths.py.

Background
----------
In amplicon sequencing (target enrichment), reads reflect the actual allele
composition of the sample. For each locus:
  - Diploid:   up to 2 alleles  (ratios expected near 0.50 / 0.50 or 1.0)
  - Tetraploid: up to 4 alleles (ratios expected near 0.25 / 0.25 / 0.25 / 0.25,
                                 or any simplex/duplex/triplex/quadruplex combination)

Dosage estimation:
  We use the fraction of reads supporting each allele peak as a proxy for
  copy number. Peaks are identified using a prominence-based peak picker on
  the read-count histogram, then dosage is assigned by comparing peak fractions
  to expected ratios (1/2 for 2x; 1/4 for 4x), with a noise floor filter.

Algorithm:
  1. Filter: keep only length values > noise floor (default 5% of total reads).
  2. Peak detection: find local maxima in the read-count histogram with a
     configurable minimum prominence (separates true alleles from stutter).
  3. Ploidy assignment: compare the number of peaks and their fractional ratios.
  4. Dosage: for each allele peak, dosage = round(fraction × ploidy).
  5. Genotype string: e.g. "120/120" (diploid homo), "114/120" (diploid het),
     "120/120/126/132" (tetraploid), "120[2]/132[2]" (tetraploid duplex).

Usage:
    python3 04_call_dosage.py \
        --dist       str_lengths/length_distributions.tsv \
        --ploidy-map ploidy.tsv \
        --outdir     calls/ \
        [--noise-floor 0.05] [--min-stutter-sep 1] [--min-reads 10]

ploidy.tsv: two-column TSV (no header) with sample name and ploidy (2 or 4).
    If a sample is not listed, it defaults to --default-ploidy (2).

Output:
    calls/allele_calls.tsv    : one row per sample-locus, allele lengths + dosages
    calls/genotypes.csv       : wide matrix (samples × loci) with genotype strings
    calls/qc_flags.tsv        : loci flagged for low coverage or ambiguous calls

Requirements: pandas, numpy, scipy
"""

import argparse
import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dist",          required=True,  help="length_distributions.tsv from script 03")
    p.add_argument("--ploidy-map",    default=None,   help="TSV: sample<TAB>ploidy (2 or 4)")
    p.add_argument("--default-ploidy",type=int, default=2, help="Ploidy for samples not in --ploidy-map")
    p.add_argument("--outdir",        default="calls")
    p.add_argument("--noise-floor",   type=float, default=0.05,
                   help="Minimum fraction of total reads for a peak to be called (default 0.05)")
    p.add_argument("--min-reads",     type=int, default=10,
                   help="Minimum total read pairs to attempt calling (default 10)")
    p.add_argument("--stutter-window",type=int, default=1,
                   help="Motif units below a major peak considered stutter (default 1)")
    p.add_argument("--dosage-tol",    type=float, default=0.15,
                   help="Tolerance around expected dosage fraction (default 0.15)")
    return p.parse_args()


# ---------------------------------------------------------------------------
def load_ploidy_map(path, default_ploidy):
    pm = {}
    if path and os.path.exists(path):
        df = pd.read_csv(path, sep="\t", header=None, names=["sample","ploidy"])
        pm = dict(zip(df["sample"], df["ploidy"].astype(int)))
        print(f"[04] Loaded ploidy for {len(pm)} samples from {path}")
    return pm, default_ploidy


# ---------------------------------------------------------------------------
def pick_peaks(lengths, counts, motif_len, noise_frac, stutter_window):
    """
    Identify allele peaks from a length:count Series.
    Returns list of (length, count, fraction) for valid peaks.
    """
    total = counts.sum()
    if total == 0:
        return []

    fracs = counts / total
    # Build dense array from min to max length
    l_min, l_max = lengths.min(), lengths.max()
    dense_len = np.arange(l_min, l_max + 1, 1)
    dense_cnt = np.zeros(len(dense_len), dtype=float)
    for l, c in zip(lengths, counts):
        dense_cnt[l - l_min] = c

    # Find peaks with minimum prominence
    prominence_threshold = noise_frac * total
    peak_idx, properties = find_peaks(
        dense_cnt,
        prominence=prominence_threshold,
        distance=motif_len  # peaks must be at least 1 motif unit apart
    )

    peaks = []
    for idx in peak_idx:
        length_val = dense_len[idx]
        cnt_val    = dense_cnt[idx]
        frac_val   = cnt_val / total

        if frac_val < noise_frac:
            continue

        peaks.append((int(length_val), int(cnt_val), round(float(frac_val), 4)))

    # Also include any length above noise floor that wasn't caught by peak finder
    # (handles cases where only 1 allele is present — flat distribution)
    if not peaks:
        for l, f in zip(lengths, fracs):
            if f >= noise_frac:
                peaks.append((int(l), int(counts[lengths == l].iloc[0]), round(float(f), 4)))

    # Remove stutter: if a peak is exactly stutter_window * motif_len below
    # a larger peak and has < 20% the reads, discard it
    peaks = sorted(peaks, key=lambda x: -x[1])  # sort by count desc
    kept = []
    for pk in peaks:
        is_stutter = False
        for major in kept:
            if (major[0] - pk[0]) in range(motif_len, (stutter_window + 1) * motif_len + 1, motif_len):
                if pk[1] < 0.20 * major[1]:
                    is_stutter = True
                    break
        if not is_stutter:
            kept.append(pk)

    return sorted(kept, key=lambda x: x[0])  # sort by length


# ---------------------------------------------------------------------------
def assign_dosage(peaks, ploidy, tol):
    """
    Given a list of (length, count, fraction) peaks and the sample ploidy,
    assign dosage (copy number) to each allele.
    Returns list of (length, dosage).
    """
    if not peaks:
        return []

    total_frac = sum(p[2] for p in peaks)
    dosages = []

    for length, cnt, frac in peaks:
        # Expected fraction = dosage / ploidy
        norm_frac = frac / total_frac if total_frac > 0 else frac
        raw_dose  = norm_frac * ploidy
        dose      = int(round(raw_dose))
        dose      = max(1, min(dose, ploidy))
        dosages.append((length, dose))

    # Ensure total dosage == ploidy; adjust largest-residual allele if needed
    total_dose = sum(d for _, d in dosages)
    if total_dose != ploidy and dosages:
        diff = ploidy - total_dose
        # Add/subtract from the allele whose raw dosage is furthest from integer
        fracs  = [p[2] / total_frac if total_frac > 0 else p[2] for p in peaks]
        raw_ds = [f * ploidy for f in fracs]
        residuals = [abs(r - round(r)) for r in raw_ds]
        adj_idx = int(np.argmax(residuals))
        dosages[adj_idx] = (dosages[adj_idx][0], max(1, dosages[adj_idx][1] + diff))

    return dosages


# ---------------------------------------------------------------------------
def format_genotype(dosages, ploidy):
    """Format dosage list into a genotype string."""
    if not dosages:
        return "./."

    # Expand alleles by dosage: e.g. [(120, 2), (132, 2)] -> "120/120/132/132"
    alleles = []
    for length, dose in sorted(dosages, key=lambda x: x[0]):
        alleles.extend([str(length)] * dose)

    if len(alleles) < ploidy:
        alleles += ["."] * (ploidy - len(alleles))

    return "/".join(alleles)


# ---------------------------------------------------------------------------
def call_sample_locus(group, ploidy, noise_floor, stutter_window, min_reads, dosage_tol):
    """Process one (sample, locus) group. Returns dict of call results."""
    motif_len   = int(group["motif_len"].iloc[0])
    total_reads = int(group["total_pairs"].iloc[0])  # same for all rows in group

    result = {
        "total_reads": total_reads,
        "n_peaks":     0,
        "alleles":     None,
        "dosages":     None,
        "genotype":    "./.",
        "flag":        "",
    }

    if total_reads < min_reads:
        result["flag"] = f"LOW_COV({total_reads})"
        return result

    lengths = group["str_length"].values
    counts  = group["read_count"].values
    # Build Series for peak picker
    len_ser = pd.Series(counts, index=lengths)
    len_ser = len_ser.sort_index()

    peaks = pick_peaks(
        len_ser.index.values,
        len_ser.values,
        motif_len,
        noise_floor,
        stutter_window
    )

    if not peaks:
        result["flag"] = "NO_PEAKS"
        return result

    # Warn if peaks > ploidy (possible contamination or paralog)
    flag = ""
    if len(peaks) > ploidy:
        flag = f"EXCESS_PEAKS({len(peaks)}>ploidy{ploidy})"
        # Keep only top-N peaks by count
        peaks = sorted(peaks, key=lambda x: -x[1])[:ploidy]
        peaks = sorted(peaks, key=lambda x: x[0])

    dosages = assign_dosage(peaks, ploidy, dosage_tol)
    genotype = format_genotype(dosages, ploidy)

    result.update({
        "n_peaks":  len(peaks),
        "alleles":  ",".join(str(p[0]) for p in peaks),
        "dosages":  ",".join(str(d[1]) for d in dosages),
        "genotype": genotype,
        "flag":     flag,
        "peak_fractions": ",".join(str(p[2]) for p in peaks),
    })
    return result


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load distributions
    dist = pd.read_csv(args.dist, sep="\t")
    print(f"[04] Loaded {len(dist)} rows from {args.dist}")

    samples = dist["sample"].unique()
    loci    = dist["locus"].unique()
    print(f"[04] {len(samples)} samples × {len(loci)} loci")

    # Load ploidy
    ploidy_map, default_ploidy = load_ploidy_map(args.ploidy_map, args.default_ploidy)

    # ---- per sample-locus call --------------------------------------------
    call_rows = []
    qc_rows   = []

    for (sample, locus), group in dist.groupby(["sample","locus"]):
        ploidy = ploidy_map.get(sample, default_ploidy)

        res = call_sample_locus(
            group,
            ploidy         = ploidy,
            noise_floor    = args.noise_floor,
            stutter_window = args.stutter_window,
            min_reads      = args.min_reads,
            dosage_tol     = args.dosage_tol,
        )

        row = {
            "sample":         sample,
            "locus":          locus,
            "ploidy":         ploidy,
            "total_reads":    res["total_reads"],
            "n_peaks":        res["n_peaks"],
            "alleles":        res.get("alleles", ""),
            "dosages":        res.get("dosages", ""),
            "genotype":       res["genotype"],
            "peak_fractions": res.get("peak_fractions", ""),
            "flag":           res["flag"],
        }
        call_rows.append(row)

        if res["flag"]:
            qc_rows.append({"sample": sample, "locus": locus, "flag": res["flag"]})

    calls_df = pd.DataFrame(call_rows)

    # ---- allele_calls.tsv -------------------------------------------------
    calls_path = os.path.join(args.outdir, "allele_calls.tsv")
    calls_df.to_csv(calls_path, sep="\t", index=False)
    print(f"[04] Wrote {calls_path}")

    # ---- genotypes.csv (wide pivot) --------------------------------------
    geno_wide = calls_df.pivot_table(
        index="sample", columns="locus", values="genotype", aggfunc="first"
    )
    geno_wide.columns.name = None
    geno_path = os.path.join(args.outdir, "genotypes.csv")
    geno_wide.to_csv(geno_path)
    print(f"[04] Wrote {geno_path}  ({geno_wide.shape[0]} samples × {geno_wide.shape[1]} loci)")

    # ---- dosage_matrix.csv -----------------------------------------------
    # For downstream popgen: rows=samples, cols=loci, value="allele[dose]" pairs
    dose_wide = calls_df.pivot_table(
        index="sample", columns="locus", values="alleles", aggfunc="first"
    )
    dose_path = os.path.join(args.outdir, "dosage_matrix.csv")
    dose_wide.to_csv(dose_path)
    print(f"[04] Wrote {dose_path}")

    # ---- qc_flags.tsv ----------------------------------------------------
    qc_path = os.path.join(args.outdir, "qc_flags.tsv")
    pd.DataFrame(qc_rows).to_csv(qc_path, sep="\t", index=False)
    print(f"[04] Wrote {qc_path}  ({len(qc_rows)} flagged sample-locus pairs)")

    # ---- summary ---------------------------------------------------------
    n_called = (calls_df["genotype"] != "./.").sum()
    n_total  = len(calls_df)
    print(f"\n[04] Call rate: {n_called}/{n_total} ({100*n_called/n_total:.1f}%)")
    if len(qc_rows):
        flag_counts = pd.DataFrame(qc_rows)["flag"].str.extract(r'^([A-Z_]+)')[0].value_counts()
        print("[04] Flag summary:")
        for flag, cnt in flag_counts.items():
            print(f"     {flag}: {cnt}")

    print(f"\n[04] Next: run  05_summarize_results.py  --calls {calls_path}")


if __name__ == "__main__":
    main()
