# Ulmus STR Amplicon Pipeline

A reference-anchored pipeline for calculating microsatellite allele lengths and
dosage from paired-end Illumina amplicon data in a mixed-ploidy system (diploid + tetraploid *Ulmus americana*).

---

## Quick note:
Hi, Marne here. I don't condone the use of AI for many things; however, I tried a lot of different microsatellite analysis programs and none of them worked with my dataset. I explained everything I needed to Claude the AI, as well as the programs I tried and why they failed. It created all of this, then I worked to troubleshoot it. It is still in progress so I will be back to update depending on how it goes.

## Why this pipeline?

| Tool tried | Why it failed here |
|---|---|
| seq2sat | STR region exceeds single read length (amplicons 94–282 bp > 150 bp read) |
| STRait Razor 3.0 | Same read-spanning problem |
| gangSTR | Size-model calibration errors on these amplicons |
| hipSTR | Stutter model requires larger flanking regions than are present |
| RepeatSeq | Obsolete dependencies, no longer installable |

**This pipeline avoids all of the above by:**
- Aligning **paired reads** to a reference with BWA-MEM2, so paired-end insert
  size bridges STR regions longer than a single read.
- Computing STR length from **insert size minus non-STR sequence**, not from
  read-level repeat counting.
- Using **pileup-based peak detection** rather than stutter models, making it
  agnostic to flanking region length.
- Handling **both 2x and 4x ploidy** explicitly with dosage assignment per locus.

---

## Requirements

### Software (must be in `$PATH`)
```
bwa-mem2 >= 2.2      (or bwa >= 0.7.17 — set USE_BWA_MEM=1 in 02_align.sh)
samtools >= 1.15
python3 >= 3.8
```

### Python packages
```bash
pip install pysam pandas numpy scipy
```

---

## Input files

| File | Description |
|---|---|
| `microsats_all_info.csv` | Microsatellite metadata CSV (provided) |
| `amplicon_reference.fa` | Reference FASTA of gene region sequences |
| `samples.txt` | One sample name per line |
| `trimmed_reads/` | Directory with `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz` |
| `ploidy.tsv` | (Optional) Two columns: `sample<TAB>ploidy` (2 or 4) |

**ploidy.tsv example:**
```
sample_001	2
sample_002	4
sample_003	4
```
Samples not listed use `--default-ploidy` (default: 2).

**To generate the amplicon_reference.fa from the microsats_all_info.csv**
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('microsats_all_info.csv').dropna(subset=['gene'])
with open('amplicon_reference.fa', 'w') as fh:
    for _, row in df.iterrows():
        fh.write(f'>{row[\"asssay_seq_name\"]}\n{row[\"assay_sequence\"]}\n')
print(f'Wrote {len(df)} amplicon sequences to amplicon_reference.fa')
"
```
---

## Usage

### Quick start (all steps)
```bash
# 1. Edit CONFIGURATION block in run_pipeline.sh
nano run_pipeline.sh

# 2. Run
bash run_pipeline.sh
```

### Step-by-step
```bash
# Step 1: Parse CSV, generate BED and config files
python3 01_prepare_loci.py \
    --csv   Ulmus_amplicon_targets.csv \
    --outdir loci/

# Step 2: Align trimmed reads
bash 02_align.sh \
    -r ulmus_reference.fa \
    -s samples.txt \
    -i trimmed_reads/ \
    -o alignments/ \
    -b loci/str_regions_padded.bed \
    -t 16

# Step 3: Extract STR lengths from insert sizes
python3 03_extract_str_lengths.py \
    --bam-dir  alignments/ \
    --config   loci/loci_config.tsv \
    --outdir   str_lengths/

# Step 4: Call alleles and dosage
python3 04_call_dosage.py \
    --dist         str_lengths/length_distributions.tsv \
    --ploidy-map   ploidy.tsv \
    --default-ploidy 2 \
    --outdir       calls/

# Step 5: Final QC and output
python3 05_summarize_results.py \
    --calls   calls/allele_calls.tsv \
    --config  loci/loci_config.tsv \
    --outdir  results/
```

---

## Output files

```
results/
├── final_alleles.csv       ← PRIMARY OUTPUT: one row per sample-locus
├── genotype_matrix.csv     ← Wide pivot (samples × loci) for popgen tools
├── locus_stats.csv         ← Per-locus call rate, allele richness, coverage
├── sample_stats.csv        ← Per-sample call rate and coverage
└── pipeline_report.txt     ← Human-readable summary

calls/
├── allele_calls.tsv        ← Detailed calls with peak fractions and flags
├── dosage_matrix.csv       ← Raw allele lengths per sample-locus
└── qc_flags.tsv            ← Flagged sample-locus pairs

str_lengths/
├── allele_lengths_raw.tsv  ← One row per spanning read pair
└── length_distributions.tsv← Aggregated length histograms

loci/
├── loci_config.tsv         ← Full locus metadata
├── loci.bed                ← STR coordinates (BED4)
├── str_regions_padded.bed  ← Padded amplicon regions for read extraction
└── primers.tsv             ← Primer sequences

alignments/
└── {sample}/
    ├── {sample}.full.bam   ← Full genome alignment
    ├── {sample}.str.bam    ← STR-region subset
    └── {sample}.str.flagstat
```

### final_alleles.csv columns

| Column | Description |
|---|---|
| sample | Sample identifier |
| locus | Microsatellite locus (gene name from CSV) |
| ploidy | Assigned ploidy (2 or 4) |
| genotype | Alleles as slash-separated lengths, e.g. `120/132` or `120/120/132/144` |
| allele_1..allele_N | Individual allele lengths (N = max ploidy) |
| dosages | Copy number per allele, e.g. `1,1` or `2,2` |
| n_peaks | Number of allele peaks detected |
| total_reads | Total spanning read pairs for this locus |
| peak_fractions | Fraction of reads per allele peak |
| flag | QC flags: LOW_COV, NO_PEAKS, EXCESS_PEAKS |

---

## Dosage calling logic

For **diploid** samples:
- 1 peak → homozygous (dosage 2/2): e.g. `120/120`
- 2 peaks at ~50:50 → heterozygous (dosage 1/1): e.g. `120/132`

For **tetraploid** samples:
- 1 peak → quadruplex (4 copies): e.g. `120/120/120/120`
- 2 peaks at ~75:25 → triplex/simplex: e.g. `120/120/120/132`
- 2 peaks at ~50:50 → duplex/duplex: e.g. `120/120/132/132`
- 3 peaks at ~50:25:25 → duplex/simplex/simplex
- 4 peaks at ~25:25:25:25 → simplex each

Dosage assignment uses the fractional read support per peak normalized to
the total and rounded to the nearest integer copy number. A tolerance of
±15% (adjustable via `--dosage-tol`) handles coverage imbalance between alleles.

---

## Tuning parameters

| Script | Parameter | Default | When to change |
|---|---|---|---|
| 01 | `--padding` | 50 bp | Increase for very long amplicons |
| 02 | `SORT_MEM` | 2G | Increase for large genomes |
| 03 | `--min-mapq` | 20 | Lower for divergent samples |
| 03 | `--min-pairs` | 5 | Increase for higher confidence |
| 03 | `--max-insert` | 1000 | Adjust to max amplicon size |
| 04 | `--noise-floor` | 0.05 | Lower to detect rare alleles |
| 04 | `--stutter-window` | 1 | Increase if stutter is severe |
| 04 | `--dosage-tol` | 0.15 | Increase for uneven coverage |
| 05 | `--min-call-rate` | 0.70 | Raise to exclude poor samples |
| 05 | `--min-locus-call-rate` | 0.50 | Raise to exclude poor loci |

---

## Downstream analysis

The `genotype_matrix.csv` (slash-separated allele lengths) is directly
compatible with:
- **polysat** (R): `read.GeneMapper()` or manual import
- **pegas** (R): `read.loci()`
- **poppr** (R): after conversion with `df2genind()`

For **Structure** / **STRUCTURE** format, each allele position needs to be
in a separate column — use the `allele_1`, `allele_2`, ... columns in
`final_alleles.csv` and convert missing data to `-9`.
