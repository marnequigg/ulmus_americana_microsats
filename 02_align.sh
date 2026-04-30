#!/usr/bin/env bash
# =============================================================================
# 02_align.sh
# =============================================================================
# Align trimmed paired-end reads to a reference genome, then sort, index, and
# produce per-sample BAM files restricted to the STR amplicon regions.
#
# This script is intentionally reference-based (not de-novo assembly) because:
#   • Many STR amplicons exceed 150 bp (single Illumina read length), so
#     individual reads cannot span the full locus. BWA-MEM2 handles paired
#     reads and allows insert-size estimation to span the gap.
#   • Reference-anchored coordinates are required for the length calculation
#     step (script 03).
#
# Tools required (must be in PATH):
#   bwa-mem2 >= 2.2        (or bwa mem as fallback; see USE_BWA_MEM flag)
#   samtools >= 1.15
#
# Usage:
#   bash 02_align.sh \
#       -r reference.fa \
#       -s sample_list.txt \
#       -i /path/to/trimmed_reads/ \
#       -o alignments/ \
#       -b loci/str_regions_padded.bed \
#       [-t 8] [-p 2x150]
#
# sample_list.txt: one sample name per line (no extensions)
#   Reads expected as: {reads_dir}/{sample}_R1.fastq.gz  {sample}_R2.fastq.gz
#   (edit READS_SUFFIX_* below if your naming differs)
#
# =============================================================================
set -euo pipefail

# ---------- defaults ---------------------------------------------------------
THREADS=8
READS_SUFFIX_R1="_R1_paired.fastq.gz"
READS_SUFFIX_R2="_R2_paired.fastq.gz"
USE_BWA_MEM=0          # set to 1 to use classic bwa mem instead of bwa-mem2
SORT_MEM="2G"

# ---------- parse args -------------------------------------------------------
usage() {
  grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,2\}//'
  exit 1
}

REF=""
SAMPLE_LIST=""
READS_DIR=""
OUT_DIR=""
BED=""

while getopts "r:s:i:o:b:t:h" opt; do
  case $opt in
    r) REF="$OPTARG" ;;
    s) SAMPLE_LIST="$OPTARG" ;;
    i) READS_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    b) BED="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$REF" || -z "$SAMPLE_LIST" || -z "$READS_DIR" || -z "$OUT_DIR" || -z "$BED" ]] && {
  echo "ERROR: -r, -s, -i, -o, and -b are all required." >&2; usage
}

# ---------- index reference if needed ----------------------------------------
echo "[02] Checking reference index..."
if [[ "$USE_BWA_MEM" -eq 1 ]]; then
  ALIGNER="bwa mem"
  [[ -f "${REF}.bwt" ]] || { echo "[02] Indexing with bwa..."; bwa index "$REF"; }
else
  ALIGNER="bwa-mem2 mem"
  [[ -f "${REF}.bwt.2bit.64" ]] || { echo "[02] Indexing with bwa-mem2..."; bwa-mem2 index "$REF"; }
fi

# samtools fai
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

mkdir -p "$OUT_DIR"

# ---------- per-sample alignment ---------------------------------------------
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
  [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue

  R1="${READS_DIR}/${SAMPLE}${READS_SUFFIX_R1}"
  R2="${READS_DIR}/${SAMPLE}${READS_SUFFIX_R2}"

  if [[ ! -f "$R1" ]]; then
    echo "[02] WARNING: R1 not found for $SAMPLE ($R1), skipping." >&2; continue
  fi
  if [[ ! -f "$R2" ]]; then
    echo "[02] WARNING: R2 not found for $SAMPLE ($R2), skipping." >&2; continue
  fi

  SAMPLE_DIR="${OUT_DIR}/${SAMPLE}"
  mkdir -p "$SAMPLE_DIR"

  FULL_BAM="${SAMPLE_DIR}/${SAMPLE}.full.bam"
  STR_BAM="${SAMPLE_DIR}/${SAMPLE}.str.bam"

  echo "[02] ---- $SAMPLE ----"

  # ---- align + sort --------------------------------------------------------
  # -R sets read group (required for downstream tools)
  # -Y  soft-clip supplementary alignments (keeps paired info intact)
  # No -F filtering here — keep all alignments; script 03 applies its own
  # quality thresholds per read.
  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}"

  echo "[02]   Aligning $SAMPLE..."
  $ALIGNER \
      -t "$THREADS" \
      -Y \
      -R "$RG" \
      "$REF" "$R1" "$R2" \
      -o "${SAMPLE_DIR}/${SAMPLE}.unsorted.bam"

  samtools sort \
      -@ "$THREADS" \
      -m "4G" \
      -o "$FULL_BAM" \
      -T "${SAMPLE_DIR}/tmp_sort" \
      "${SAMPLE_DIR}/${SAMPLE}.unsorted.bam"

  rm "${SAMPLE_DIR}/${SAMPLE}.unsorted.bam"

  samtools index -@ "$THREADS" "$FULL_BAM"

  # ---- subset to STR amplicon regions only ----------------------------------
  # -L restricts to BED regions; -F 2308 removes unmapped/secondary/supplementary
  # -q 20 minimum mapping quality
  echo "[02]   Subsetting to STR regions..."
  samtools view \
      -@ "$THREADS" \
      -b \
      -L "$BED" \
      -F 2308 \
      -q 20 \
      -o "$STR_BAM" \
      "$FULL_BAM"

  samtools index -@ "$THREADS" "$STR_BAM"

  # ---- flagstat for QC ------------------------------------------------------
  samtools flagstat "$STR_BAM" > "${SAMPLE_DIR}/${SAMPLE}.str.flagstat"
  echo "[02]   Done: $STR_BAM"

done < "$SAMPLE_LIST"

echo ""
echo "[02] All samples aligned. BAMs in: $OUT_DIR"
echo "[02] Next: run  03_extract_str_lengths.py  --bam-dir $OUT_DIR  --config loci/loci_config.tsv"
