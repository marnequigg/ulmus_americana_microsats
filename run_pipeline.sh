#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  —  Ulmus STR amplicon pipeline (master runner)
# =============================================================================
# Executes all five scripts in order. Edit the CONFIGURATION block below
# before running.
#
# Usage:
#   bash run_pipeline.sh
#
# Logs go to: logs/pipeline_YYYYMMDD_HHMMSS.log
# =============================================================================
set -euo pipefail

# =============================================================================
# CONFIGURATION — edit these paths before running
# =============================================================================

# Path to the Ulmus amplicon microsatellite CSV
CSV="/home/imh4101/05b.microsats_claude/microsats_all_info.csv"

# Reference genome FASTA (must be indexed with bwa-mem2; script 02 will index if absent)
REFERENCE="/home/imh4101/05b.microsats_claude/amplicon_reference.fa"

# One sample name per line (no file extensions)
SAMPLE_LIST="/home/imh4101/05b.microsats_claude/samples.txt"

# Directory containing trimmed paired-end reads
# Reads expected: {READS_DIR}/{sample}_R1.fastq.gz  {sample}_R2.fastq.gz
READS_DIR="/home/imh4101/05b.microsats_claude/01.input_trimmed_reads/"

# Two-column TSV (sample<TAB>ploidy) — omit samples for default ploidy below
PLOIDY_MAP="/home/imh4101/05b.microsats_claude/ploidy.tsv"
DEFAULT_PLOIDY=2          # used for samples not in PLOIDY_MAP

# CPUs for alignment
THREADS=30

# Output subdirectories (created automatically)
LOCI_DIR="loci"
ALIGN_DIR="alignments"
LENGTHS_DIR="str_lengths"
CALLS_DIR="calls"
RESULTS_DIR="results"
LOG_DIR="logs"

# =============================================================================
# END CONFIGURATION
# =============================================================================

mkdir -p "$LOG_DIR"
LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo " Ulmus STR Amplicon Pipeline"
echo " Started: $(date)"
echo "============================================================"
echo ""

# ---- check dependencies ----------------------------------------------------
echo "[run] Checking dependencies..."
MISSING=0
for cmd in python3 bwa-mem2 samtools; do
  if ! command -v "$cmd" &>/dev/null; then
    echo "  MISSING: $cmd"
    MISSING=1
  else
    echo "  OK: $cmd ($(command -v $cmd))"
  fi
done

python3 -c "import pysam, pandas, numpy, scipy" 2>/dev/null && \
  echo "  OK: python packages (pysam, pandas, numpy, scipy)" || \
  { echo "  MISSING python packages. Install with:"; \
    echo "    pip install pysam pandas numpy scipy"; MISSING=1; }

if [[ "$MISSING" -eq 1 ]]; then
  echo "[run] ERROR: Missing dependencies. Aborting."
  exit 1
fi
echo ""

# ---- Step 1: Prepare loci --------------------------------------------------
echo "------------------------------------------------------------"
echo "[run] STEP 1: Prepare loci"
echo "------------------------------------------------------------"
python3 /home/imh4101/05b.microsats_claude/scripts/01_prepare_loci.py \
    --csv    "$CSV" \
    --outdir "$LOCI_DIR"
echo ""

# ---- Step 2: Align ---------------------------------------------------------
echo "------------------------------------------------------------"
echo "[run] STEP 2: Align reads to reference"
echo "------------------------------------------------------------"
bash /home/imh4101/05b.microsats_claude/scripts/02_align.sh \
    -r "$REFERENCE" \
    -s "$SAMPLE_LIST" \
    -i "$READS_DIR" \
    -o "$ALIGN_DIR" \
    -b "${LOCI_DIR}/str_regions_padded.bed" \
    -t "$THREADS"
echo ""

# ---- Step 3: Extract STR lengths -------------------------------------------
echo "------------------------------------------------------------"
echo "[run] STEP 3: Extract STR lengths from BAMs"
echo "------------------------------------------------------------"
python3 /home/imh4101/05b.microsats_claude/scripts/03_extract_str_lengths.py \
    --bam-dir  "$ALIGN_DIR" \
    --config   "${LOCI_DIR}/loci_config.tsv" \
    --outdir   "$LENGTHS_DIR"
echo ""

# ---- Step 4: Call alleles and dosage ---------------------------------------
echo "------------------------------------------------------------"
echo "[run] STEP 4: Call alleles and dosage"
echo "------------------------------------------------------------"
PLOIDY_ARG=""
[[ -f "$PLOIDY_MAP" ]] && PLOIDY_ARG="--ploidy-map $PLOIDY_MAP"

python3 /home/imh4101/05b.microsats_claude/scripts/04_call_dosage.py \
    --dist          "${LENGTHS_DIR}/length_distributions.tsv" \
    $PLOIDY_ARG \
    --default-ploidy "$DEFAULT_PLOIDY" \
    --outdir         "$CALLS_DIR"
echo ""

# ---- Step 5: Summarize results ---------------------------------------------
echo "------------------------------------------------------------"
echo "[run] STEP 5: Summarize results"
echo "------------------------------------------------------------"
python3 /home/imh4101/05b.microsats_claude/scripts/05_summarize_results.py \
    --calls   "${CALLS_DIR}/allele_calls.tsv" \
    --config  "${LOCI_DIR}/loci_config.tsv" \
    --outdir  "$RESULTS_DIR"
echo ""

echo "============================================================"
echo " Pipeline complete: $(date)"
echo " Primary output: ${RESULTS_DIR}/final_alleles.csv"
echo " Full log:       $LOG"
echo "============================================================"
