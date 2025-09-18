#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "Dataset\tLanguage\tRuntime\tN50"
echo -e "-------------------------------------"

for dataset in data{1..4}; do

    # Run Python Script
    start=$(date +%s)
    result=$(python3 "${DIR}/main.py" "${DIR}/${dataset}")
    end=$(date +%s)
    runtime=$((end - start))
    runtime_fmt=$(date -ud "@$runtime" +'%H:%M:%S')

    n50=$(echo "$result"| awk '{print $2}')
    echo -e "${dataset}\tPython\t${runtime_fmt}\t${n50}"    start=$(date +%s)

    # Run Codon Script
    start=$(date +%s)
    result=$(codon run -release "main_codon.py" "/${dataset}")
    end=$(date +%s)
    runtime=$((end - start))
    runtime_fmt=$(date -ud "@$runtime" +'%H:%M:%S')

    n50=$(echo "$result"| awk '{print $2}')
    echo -e "${dataset}\tCodon\t${runtime_fmt}\t${n50}"

done