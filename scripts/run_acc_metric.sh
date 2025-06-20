#!/usr/bin/env bash
set -euo pipefail

# Path to the Python evaluation script
EVAL_SCRIPT="/mnt/data/luoling/benchmark_docker/131Ab_hcd_denovo/acc_metric.py"

# Root directory containing enzyme_results subfolders
ROOT_DIR="/mnt/data/luoling/benchmark_docker/131Ab_hcd_denovo/enzyme_results"

# Iterate over every .csv file except any acc_metric.csv or all_results.csv
find "$ROOT_DIR" -type f -name '*.csv' \
     ! -name 'acc_metric.csv' \
     ! -name 'all_results.csv' | while read -r INPUT_CSV; do

  # Determine the directory containing this CSV
  DIR=$(dirname "$INPUT_CSV")

  # Define the output metrics file
  OUTPUT_CSV="$DIR/acc_metric.csv"

  echo "▶ Running evaluation on: $INPUT_CSV"
  echo "  → Writing metrics to: $OUTPUT_CSV"

  # Invoke the Python script
  python3 "$EVAL_SCRIPT" \
    --input "$INPUT_CSV" \
    --output "$OUTPUT_CSV"

done

echo "✅ All evaluations complete."