#!/usr/bin/env bash
# run_metrics.sh
# Usage: ./run_metrics.sh <mgf_in_path> <csv_in_path> <csv_out_path>

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <mgf_in> <csv_in> <csv_out>"
  exit 1
fi

MGF_IN="$1"
CSV_IN="$2"
CSV_OUT="$3"

# Activate virtualenv if needed:
# source /path/to/venv/bin/activate

echo "🔹 Input MGF:  $MGF_IN"
echo "🔹 Input CSV:  $CSV_IN"
echo "🔹 Output CSV: $CSV_OUT"
echo

python3 compute_metrics.py \
  --mgf-in "$MGF_IN" \
  --csv-in "$CSV_IN" \
  --csv-out "$CSV_OUT"

echo
echo "✅ Done: results written to $CSV_OUT"
