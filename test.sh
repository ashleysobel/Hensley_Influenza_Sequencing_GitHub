#!/usr/bin/env bash
set -Eeuo pipefail

# paths relative to repo root (mounted into the container)
OUTDIR="workTest"
FASTQ_DIR="test/Dino_Fastq"
RUNS="BZPB3P"

# Clean + recreate output directory
if [[ -d "$OUTDIR" ]]; then
  echo "Cleaning existing $OUTDIR/"
  rm -rf "$OUTDIR"
fi
mkdir -p "$OUTDIR"

# Use the IRMA installed in the image (on PATH)
IRMA_BIN="$(command -v IRMA)"
if [[ -z "$IRMA_BIN" ]]; then
  echo "IRMA not found on PATH"; exit 1
fi

# Run your pipeline inside the container's environment
Rscript main.R \
  -o "$OUTDIR" \
  -f "$FASTQ_DIR" \
  --runs "$RUNS" \
  --min_cov 5 
  
cat "$OUTDIR/BZPB3P_coverage_summary.csv" || echo "output csv not found"


# Helper: compare two CSVs ignoring quotes and line order
compare_csv() {
  local produced="$1"
  local expected="$2"
  if diff <(sed 's/"//g' "$produced" | sort) \
          <(sed 's/"//g' "$expected" | sort) >/dev/null; then
    echo "✅ $(basename "$produced") matches expected (ignoring quotes)"
  else
    echo "❌ $(basename "$produced") differs from expected"
    diff -u <(sed 's/"//g' "$produced" | sort) \
           <(sed 's/"//g' "$expected" | sort) || true
    exit 1
  fi
}

# Run comparisons (adjust file names)
compare_csv "$OUTDIR/BZPB3P_coverage_summary.csv" "test/BZPB3P_coverage_summary.csv"


echo "✅ All tests passed"



